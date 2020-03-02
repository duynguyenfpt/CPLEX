import cplex

model = cplex.Cplex()

variables = ['x1','x2']

model.objective.set_sense(model.objective.sense.minimize)
# model.variables.add(names=variables,types=[model.variables.type.binary,model.variables.type.binary])
model.variables.add(names=variables)

quad_variables_coefficient = [[[0,1],[1,0]],[[0,1],[0,1]]]
# for i in range(len(variables)):
#     for i2 in range(len(variables)):
#         quad_variables.append((i,i2,1))

# print(quad_variables_coefficient)

# model.objective.set_quadratic_coefficients(quad_variables_coefficient)
model.objective.set_quadratic(quad_variables_coefficient)
# model.objective.set_quadratic_coefficients([("x2","x2",1),("x1","x1",1)])
print(model.objective.get_quadratic())
# model.objective.set_linear([(0,1),(1,1)])
# model.quadratic_constraints.add(name="quad_constraint",quad_expr=cplex.SparseTriple(ind1=[i for i in variables]
#                                                                                     ,ind2=[i for i in variables]
#                                                                                     ,val=[1 for i in range(len(variables))]),
#                                 rhs=18,
#                                 sense='L')

model.linear_constraints.add(lin_expr=[[["x1","x2"],[1,1]]],senses=["E"],rhs=[11])
model.solve()
print(model.solution.get_status())
print(model.solution.get_objective_value())

