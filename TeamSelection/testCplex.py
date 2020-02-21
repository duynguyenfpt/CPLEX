import cplex

model = cplex.Cplex()

variables = ['x1','x2']

model.objective.set_sense(model.objective.sense.minimize)
model.variables.add(names=variables)

quad_variables = [[[0, 1], [4.0, 6]],
                 [[0, 1], [6, 9]]]
# for i in range(len(variables)):
#     for i2 in range(len(variables)):
#         quad_variables.append((i,i2,1))

print(quad_variables)

model.objective.set_quadratic(quad_variables)
print(model.objective.get_quadratic())
# model.quadratic_constraints.add(name="quad_constraint",quad_expr=cplex.SparseTriple(ind1=[i for i in variables]
#                                                                                     ,ind2=[i for i in variables]
#                                                                                     ,val=[1 for i in range(len(variables))]),
#                                 rhs=18,
#                                 sense='L')
model.solve()

if (model.solution.status.infeasible):
    print("Infeasible")
else:
    print("Solution for problems using MDSB model is: ")
    print(model.solution.get_values())

