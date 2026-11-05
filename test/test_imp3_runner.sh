#!/bin/bash
# Tests of the imp3 runner without Snakemake or Slurm: stub snakemake / sbatch / nvidia-smi record their calls.
# bash test/test_imp3_runner.sh   (needs python3 with PyYAML)
set -euo pipefail
IMP3=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d); H1=$(mktemp -d); H2=$(mktemp -d); trap 'rm -rf "$T" "$H1" "$H2"' EXIT
mkdir -p "$T/bin"
cat > "$T/bin/snakemake" <<'S'
#!/bin/bash
echo "SNAKEMAKE $*" >> "$CALLS"
cfg=""; prev=""; for a in "$@"; do [ "$prev" = "--configfile" ] && cfg=$a; prev=$a; done
case " $* " in *" --unlock "*|*" --dry-run "*) exit 0;; esac
out=$(python3 -c "import yaml;print(yaml.safe_load(open('$cfg'))['outputdir'])")
mkdir -p "$out/status" && touch "$out/status/preprocessing.done" "$out/status/all.done"
S
printf '#!/bin/bash\necho "SBATCH $*" >> "$CALLS"; echo "Submitted batch job 42"\n' > "$T/bin/sbatch"
printf '#!/bin/bash\necho "GPU 0: test"\n' > "$T/bin/nvidia-smi"
chmod +x "$T/bin/"*
export PATH="$T/bin:$PATH" CALLS="$T/calls.txt"
# copy IMP3 (the container markers are written next to the Snakefile)
mkdir -p "$T/imp3" && cp -r "$IMP3/imp3" "$IMP3/profiles" "$IMP3/runscripts" "$T/imp3/" && mkdir -p "$T/imp3/workflow/envs"
cp "$IMP3"/workflow/envs/*.yaml "$T/imp3/workflow/envs/"
I="$T/imp3/imp3"
ok() { echo "ok    $1"; }
fail() { echo "FAIL  $1"; exit 1; }

# binny_mantis.cfg written by binny's setup in another checkout: paths moved to this one
B="$T/imp3/submodules/binny"; mkdir -p "$B/config" "$B/database/hmms/checkm_tf" "$B/database/hmms/checkm_pf"
touch "$B/database/hmms/checkm_tf/checkm_filtered_tf.hmm" "$B/database/hmms/checkm_pf/checkm_filtered_pf.hmm"
printf 'nog_ref_folder=NA\ncustom_ref=/old/imp3/submodules/binny/database/hmms/checkm_tf/checkm_filtered_tf.hmm\ncustom_ref=/old/imp3/submodules/binny/database/hmms/checkm_pf/checkm_filtered_pf.hmm\n' > "$B/config/binny_mantis.cfg"

# slurm project, two samples
"$I" init "$T/p" --site generic --db /data/db >/dev/null 2>&1
for s in S1 S2; do mkdir -p "$T/p/samples/$s"; touch "$T/p/samples/$s/${s}_mg_R1.fq.gz" "$T/p/samples/$s/${s}_mg_R2.fq.gz"; done
touch "$T/p/samples/S2/S2_mt_R1.fq.gz" "$T/p/samples/S2/S2_mt_R2.fq.gz" "$T/p/samples/S2/readme.txt"
"$I" run -p "$T/p" -j 2 >/dev/null 2>&1
[ -f "$T/p/imp/S1/status/all.done" ] && [ -f "$T/p/imp/S2/status/all.done" ] && ok "slurm: both samples run" || fail "slurm run"
python3 - "$T/p/config/imp_config_S2.yaml" <<'P' && ok "config: reads, outputdir, db" || fail "config"
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
assert c["raws"]["Metagenomics"].endswith("S2_mg_R1.fq.gz " + c["raws"]["Metagenomics"].split()[1].rsplit("/",1)[0] + "/S2_mg_R2.fq.gz")
assert len(c["raws"]["Metatranscriptomics"].split()) == 2 and "readme" not in str(c["raws"])
assert c["outputdir"].endswith("/p/imp/S2") and c["db_path"] == "/data/db" and c["sample"] == "S2"
P
grep -q -- "--profile $T/p/profile --conda-prefix" "$CALLS" && ok "slurm: profile + conda prefix" || fail "slurm cmd"
grep -q "^custom_ref=$B/database/hmms/checkm_tf/checkm_filtered_tf.hmm$" "$B/config/binny_mantis.cfg" \
  && grep -q "/old/imp3" "$B/config/binny_mantis.cfg.orig" && ok "binny_mantis.cfg: paths of this checkout" || fail "binny cfg"
"$I" run -p "$T/p" 2>&1 | grep -q "0 to run" && ok "finished samples skipped" || fail "skip"
"$I" status -p "$T/p" | grep -q "^S1 *done" && ok "status" || fail "status"
rm -rf "$T/p/imp/S1"; "$I" run -p "$T/p" --submit >/dev/null 2>&1
grep -q "SBATCH --job-name imp3_S1 .*--wrap .*imp3 run -p $T/p --samples S1 --foreground" "$CALLS" && ok "submit master" || fail "submit"
"$I" status -p "$T/p" | grep -q "^S1 *submitted" && ok "status: submitted until the master starts" || fail "status submitted"
# max_threads: cores_max of the sample config capped
"$I" run -p "$T/p" --samples S2 --force --max-threads 16 >/dev/null 2>&1
python3 - "$T/p/config/imp_config_S2.yaml" <<'P' && ok "--max-threads caps cores_max" || fail "max threads"
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))["cluster_specs"]
assert all(c[k]["cores_max"] <= 16 and c[k]["big_mem_cores_max"] <= 16 for k in ("main_cluster", "high_mem_cluster"))
P
# master as a Slurm job on another cluster: found in the user's squeue list
mkdir -p "$T/p/imp/S1"; rm -f "$T/p/imp_run_logs/imp_S1.submitted"; printf "== x\n" > "$T/p/imp_run_logs/imp_S1_098.log"
printf '#!/bin/bash\necho "CLUSTER: iris"; echo "6283057 RUNNING"\n' > "$T/bin/squeue"; chmod +x "$T/bin/squeue"
echo "otherhost 1 6283057 imp_S1_098.log" > "$T/p/imp_run_logs/imp_S1.running"
"$I" status -p "$T/p" | grep -q "^S1 *running " && ok "status: master job found with squeue -M all" || fail "status squeue"
printf '#!/bin/bash\necho "CLUSTER: iris"\n' > "$T/bin/squeue"
"$I" status -p "$T/p" | grep -q "^S1 *stopped" && ok "status: master job gone -> stopped" || fail "status squeue gone"
rm -f "$T/bin/squeue" "$T/p/imp_run_logs/imp_S1.running" "$T/p/imp_run_logs/imp_S1_098.log"
# running master (marker with host + pid), errors of keep-going, stale marker
rm -f "$T/p/imp_run_logs/imp_S1.submitted"; mkdir -p "$T/p/imp/S1"
printf "== x\nError in rule trimming:\n" > "$T/p/imp_run_logs/imp_S1_099.log"
sleep 30 & SP=$!
echo "$(python3 -c 'import socket;print(socket.gethostname())') $SP - imp_S1_099.log" > "$T/p/imp_run_logs/imp_S1.running"
"$I" status -p "$T/p" | grep -q "^S1 *running+err" && ok "status: running (with errors, keep-going)" || fail "status running"
kill $SP; wait $SP 2>/dev/null || true; while kill -0 $SP 2>/dev/null; do sleep 0.2; done
"$I" status -p "$T/p" | grep -q "^S1 *failed" && ok "status: stale marker -> failed" || fail "status stale"
rm -f "$T/p/imp_run_logs/imp_S1.running" "$T/p/imp_run_logs/imp_S1_099.log"

"$I" init "$T/u" --site ulhpc --container /img/imp3.sif --db /data/db >/dev/null 2>&1
grep -q '^tmp_root: /scratch/users/$USER/imp3_tmp' "$T/u/imp3.yaml" && grep -q "container_mounts:" "$T/u/imp3.yaml" \
  && ok "init ulhpc container: tmp_root on scratch" || fail "init ulhpc"

# local project with container
"$I" init "$T/l" --executor local --container /img/imp3.sif --db /data/db >/dev/null 2>&1
mkdir -p "$T/l/samples/A"; touch "$T/l/samples/A/mg.r1.fq.gz" "$T/l/samples/A/mg.r2.fq.gz"
"$I" run -p "$T/l" >/dev/null 2>&1
tail -1 "$CALLS" | grep -q -- "--use-conda --use-singularity --singularity-args --nv --bind /data/db:/data/db --bind $T" && ok "container args" || fail "container args"
tail -1 "$CALLS" | grep -q -- "--cores [0-9]* --resources gpu=1" && ok "local cores + gpu" || fail "local"
out=$(SLURM_CPUS_ON_NODE=1 SLURM_MEM_PER_CPU=2048 "$I" run -p "$T/l" --force -n 2>&1 || true)
grep -q "1 cores, 2 GB for IMP3" <<< "$out" && ok "local: resources of the Slurm allocation" || fail "local slurm allocation"
which() { echo "exported which"; }; export -f which
"$I" run -p "$T/l" --force -n 2>&1 | grep -q "unset -f which && snakemake" && ok "container: exported which function unset" \
  || fail "unset which"
unset -f which
# reads through a chain of links (project -> H1 -> H2): every hop's folder is bound
mkdir -p "$T/l/samples/B" "$H1/hop" "$H2/raw"; for r in r1 r2; do touch "$H2/raw/B_$r.fastq.gz"
  ln -s "$H2/raw/B_$r.fastq.gz" "$H1/hop/mg.$r.fq.gz"; ln -s "$H1/hop/mg.$r.fq.gz" "$T/l/samples/B/mg.$r.fq.gz"; done
"$I" run -p "$T/l" --samples B >/dev/null 2>&1
tail -1 "$CALLS" | grep -q -- "--bind $H1/hop:$H1/hop" && tail -1 "$CALLS" | grep -q -- "--bind $H2/raw:$H2/raw" \
  && ok "binds along a symlink chain" || fail "symlink chain binds"
# mount points of the image: folders below them are bound as the mount point
python3 - "$T/l/imp3.yaml" "$H1" <<'P'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1])); c["container_mounts"] = [sys.argv[2]]; yaml.safe_dump(c, open(sys.argv[1], "w"))
P
"$I" run -p "$T/l" --samples B --force >/dev/null 2>&1
tail -1 "$CALLS" | grep -q -- "--bind $H1:$H1" && ! tail -1 "$CALLS" | grep -q -- "--bind $H1/hop" \
  && ok "container_mounts: bound as the mount point" || fail "container_mounts"
# bind check in the image (stub apptainer): a folder missing in the container stops the sample before Snakemake
cat > "$T/bin/apptainer" <<'S'
#!/bin/bash
[ -n "${STUB_MISSING:-}" ] && echo "MISSING $STUB_MISSING"
[ -n "${STUB_PRIVATE_TMP:-}" ] && for x in "$@"; do case $x in *.imp3_check_*) echo "MISSING $x";; esac; done
echo IMP3_BINDS_CHECKED
S
chmod +x "$T/bin/apptainer"; n=$(wc -l < "$CALLS")
out=$(STUB_MISSING="$H2/raw" "$I" run -p "$T/l" --samples B --force 2>&1 || true)
grep -q "not visible in the container: $H2/raw" <<< "$out" && [ "$(wc -l < "$CALLS")" = "$n" ] \
  && ok "bind check stops the sample" || fail "bind check"
"$I" run -p "$T/l" --samples B --force >/dev/null 2>&1 && ok "bind check passes" || fail "bind check pass"
out=$(STUB_PRIVATE_TMP=1 "$I" run -p "$T/l" --samples B --force 2>&1 || true)
grep -q "temporary folder /tmp is not shared with the container" <<< "$out" && ! ls -d /tmp/.imp3_check_* 2>/dev/null \
  && ok "private /dev/shm in the container detected" || fail "tmp probe"
rm "$T/bin/apptainer"
[ -f "$T/imp3/submodules/comebin.ready" ] && [ -f "$T/imp3/workflow/envs/binny.ready" ] && ok "setup markers" || fail "markers"
python3 - "$T/l/config/imp_config_A.yaml" <<'P' && ok "local cluster_specs valid for getThreads" || fail "cluster_specs"
import sys, yaml, math
c = yaml.safe_load(open(sys.argv[1]))
assert c["container"] == "/img/imp3.sif" and c["tmp_root"] == "/tmp"
for k in ("main_cluster", "high_mem_cluster"):
    s = c["cluster_specs"][k]
    assert math.floor(s["big_mem_total"] / s["big_mem_per_core"]) == s["big_mem_cores"]
    assert s["cores_max"] <= math.floor(s["mem_total"] / s["mem_per_core"]) and s["big_mem_cores_max"] <= s["big_mem_cores"]
P

# containerize: environments found by Snakemake's hash under a prefix
P="$T/envs"; mkdir -p "$P"
python3 - "$P" "$T/imp3/workflow/envs/IMP_mapping.yaml" <<'P'
import hashlib, os, sys
p, f = sys.argv[1], sys.argv[2]
m = hashlib.md5(); m.update(os.path.realpath(p).encode()); m.update(open(f, "rb").read())
os.makedirs(os.path.join(p, m.hexdigest()[:8], "conda-meta"))
P
"$I" containerize --envs "$P" -o "$T/c" >/dev/null 2>&1
H=$(python3 -c "import hashlib;print(hashlib.md5(open('$T/imp3/workflow/envs/IMP_mapping.yaml','rb').read()).hexdigest())")
grep -q "IMP_mapping.yaml.*/conda-envs/$H" "$T/c/envs.tsv" && bash -n "$T/c/pack_envs.sh" && bash -n "$T/c/build.sh" \
    && grep -q "tar -C / -cf envs/$H.tar.part" "$T/c/pack_envs.sh" && ok "containerize" || fail "containerize"
"$I" containerize --envs "$P" -o "$T/c2" --mounts /work/projects,/scratch/users >/dev/null 2>&1
grep -q "mkdir -p /conda-envs /work/projects /scratch/users" "$T/c2/imp3.def" && ok "containerize --mounts" \
    || fail "containerize --mounts"
! "$I" containerize --envs "$P" -o "$T/c3" --mounts "$T" >/dev/null 2>&1 && ok "--mounts over the environments refused" \
    || fail "--mounts check"
# tools: programs of the rules, checked in the container (stub apptainer runs the check on this machine)
"$I" tools --imp3-dir "$IMP3" --list | grep -q "^IMP_assembly_merge.yaml.*parallel" && ok "tools --list" || fail "tools --list"
cat > "$T/bin/apptainer" <<'S'
#!/bin/bash
while [ "$1" != "exec" ]; do shift; done; shift 2; exec "$@"
S
chmod +x "$T/bin/apptainer"
out=$("$I" tools --imp3-dir "$IMP3" --container /img/imp3.sif 2>&1 || true)
grep -q "environments not in the image" <<< "$out" && grep -q "program calls checked in /img/imp3.sif" <<< "$out" \
  && ok "tools check in the container" || fail "tools check"
rm "$T/bin/apptainer"
echo "all tests passed"
