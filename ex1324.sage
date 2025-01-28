# Loading the Hypermatrix Algebra Package.
load('Hypermatrix_Algebra_Package_code.sage')

# Initialization of the number of vertices.
sz=Integer(5)

# Initialization of the list of vertex variables.
X=var_list('x',sz)

# Initialization of the functional directed graph as a list of ordered pairs.
Tf=[(Integer(0),Integer(0))]+[(i,i-Integer(1)) for i in rg(1,sz)]

# Initialization of the vertex Vandermonde factor
Vv=prod(X[v]-X[u] for u in rg(sz) for v in rg(sz) if u<v)

# Initialization of the induced edge label Vandermonde factor
Ve=prod((X[Tf[v][1]]-X[v])^2-(X[Tf[u][1]]-X[u])^2 for u in rg(sz) for v in rg(sz) if u<v)

# Computing Pf and its remainder
Pf0=Vv*Ve; cPf0=remainder_via_lagrange_interpolation(Pf0, rg(sz), X)