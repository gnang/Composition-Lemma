# Loading the Hypermatrix Algebra Package.
load('Hypermatrix_Algebra_Package_code.sage')

# Initialization of the number of vertices
sz=Integer(5)

# Initialization of the list of vertex variables
X=var_list('x',sz)

# Initialization of the functional directed graph as a list of ordered pairs
Tf=[(Integer(0),Integer(0))]+[(i,i-Integer(1)) for i in rg(1,sz)]

# Initialization of the polynomial prior to the complementary labeling involution transformation
Vv0=prod(X[v]-X[u] for u in rg(sz) for v in rg(sz) if u<v)
Ve0=prod((X[Tf[v][1]]-X[v])^2-(X[Tf[u][1]]-X[u])^2 for u in rg(sz) for v in rg(sz) if u<v)

# Computing Pf and its remainder
Pf0=Vv0*Ve0; cPf0=remainder_via_lagrange_interpolation(Pf0, rg(sz), X)

# Initialization of the polynomial after the complementary labeling involution transformation
Vv1=prod(X[v]-X[u] for u in rg(sz) for v in rg(sz) if u<v)
Ve1=prod((X[Tf[v][1]]-X[v])^2-(X[Tf[u][1]]-X[u])^2 for u in rg(sz) for v in rg(sz) if u<v)

# Computing the image of Pf and its remainder
Pf1=Vv1*Ve1; cPf1=remainder_via_lagrange_interpolation(Pf1, rg(sz), X)
