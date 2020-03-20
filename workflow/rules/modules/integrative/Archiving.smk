# include rules to archive results

arch_ctrl = ['status/workflow.done']
if config['archiving']['keep_all']:
    include:
        '../../Archiving/archiving.smk'
    arch_ctrl.append()
elif config['archiving']['reduce']:

elif config['archiving']['maxReduce_addDB']:

