

from bbpower_config.generate_config import generate_bbcompsep_dict
from sacc_object.generate_sacc import write_tracer

bool_array = [True, False]
for cmb in bool_array:
    for dust in bool_array:
        for sync in bool_array:
            for decorr in bool_array:

                if (dust == True & sync == True):
                    for cross in bool_array:
                        generate_bbcompsep_dict(cmb, dust, sync, decorr, cross)
                else:
                    generate_bbcompsep_dict(cmb, dust, sync, decorr, cross = False)

write_tracer()
