'''
theory cell
'''
import numpy as np
from QEfgs.theory.theory_cell import get_EBlm_plm_theo #,simulate_cmbalms_theory

get_EBlm_plm_theo()

# # didnt make any difference to specify seed
# a = simulate_cmbalms_theory()

# b = simulate_cmbalms_theory(10)
# c = simulate_cmbalms_theory(1)

# print(np.all(np.isclose(a,b)))
# print(np.all(np.isclose(a,c)))
