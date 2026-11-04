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
"$I" run -p "$T/p" 2>&1 | grep -q "0 to run" && ok "finished samples skipped" || fail "skip"
"$I" status -p "$T/p" | grep -q "^S1 *done" && ok "status" || fail "status"
rm -rf "$T/p/imp/S1"; "$I" run -p "$T/p" --submit >/dev/null 2>&1
grep -q "SBATCH --job-name imp3_S1 .*--wrap .*imp3 run -p $T/p --samples S1 --foreground" "$CALLS" && ok "submit master" || fail "submit"

# local project with container
"$I" init "$T/l" --executor local --container /img/imp3.sif --db /data/db >/dev/null 2>&1
grep -q "^tmp_root: /tmp" "$T/l/imp3.yaml" && sed -i "s|^tmp_root: /tmp|tmp_root: /dev/shm|" "$T/l/imp3.yaml"   # test folders are in /tmp
mkdir -p "$T/l/samples/A"; touch "$T/l/samples/A/mg.r1.fq.gz" "$T/l/samples/A/mg.r2.fq.gz"
"$I" run -p "$T/l" >/dev/null 2>&1
tail -1 "$CALLS" | grep -q -- "--use-conda --use-singularity --singularity-args --nv --bind /dev/shm:/dev/shm --bind /data/db:/data/db --bind $T" && ok "container args" || fail "container args"
tail -1 "$CALLS" | grep -q -- "--cores [0-9]* --resources gpu=1" && ok "local cores + gpu" || fail "local"
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
printf '#!/bin/bash\n[ -n "${STUB_MISSING:-}" ] && echo "MISSING $STUB_MISSING"; echo IMP3_BINDS_CHECKED\n' > "$T/bin/apptainer"
chmod +x "$T/bin/apptainer"; n=$(wc -l < "$CALLS")
out=$(STUB_MISSING="$H2/raw" "$I" run -p "$T/l" --samples B --force 2>&1 || true)
grep -q "not visible in the container: $H2/raw" <<< "$out" && [ "$(wc -l < "$CALLS")" = "$n" ] \
  && ok "bind check stops the sample" || fail "bind check"
"$I" run -p "$T/l" --samples B --force >/dev/null 2>&1 && ok "bind check passes" || fail "bind check pass"
rm "$T/bin/apptainer"
[ -f "$T/imp3/submodules/comebin.ready" ] && [ -f "$T/imp3/workflow/envs/binny.ready" ] && ok "setup markers" || fail "markers"
python3 - "$T/l/config/imp_config_A.yaml" <<'P' && ok "local cluster_specs valid for getThreads" || fail "cluster_specs"
import sys, yaml, math
c = yaml.safe_load(open(sys.argv[1]))
assert c["container"] == "/img/imp3.sif" and c["tmp_root"] == "/dev/shm"
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
echo "all tests passed"
