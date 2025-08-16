using GadgetSearch
using Gurobi

const env = Gurobi.Env()

truth_table = GadgetSearch.generic_rule(110, (3, 1))

generate_full_grid_udg(Triangular(), 4, 4)

results, failed = GadgetSearch.search_by_truth_tables(
           GraphLoader("udg.g6"), 
           [truth_table];
           optimizer =Gurobi.Optimizer, env=env, pin_candidates=[[1,2,3,4],[2,3,4,1]], allow_defect=true,  objective=x->sum(x), max_result_num=1
       )







