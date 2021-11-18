##################################################
# Mantis
##################################################
# https://git-r3lab.uni.lu/pedro.queiros/mantis
# NOTE: expecting Mantis to be included as a git submodule

# TODO

# rule mantis_config:
#     output:
#         # config file
#     log: # if needed
#         # ???
#     message:
#         "Mantis: create config file for the setup"
#     shell/run/script: # TODO
#         # TODO

# rule mantis_setup_dbs:
#     input:
#         # config file
#     output:
#         # ???
#     log:
#         # ???
#     threads: 1 # can remove if only one 1 can be used
#     params:
#         mantis_folder=os.path.join(MOD_DIR, "mantis")
#     message:
#         "Mantis: setup ???"
#     conda:
#         os.path.join(MOD_DIR, "mantis", "mantis_env.yml")
#     shell:
#         "python {params.mantis_folder} setup_databases -mc {input} &> {log}"
