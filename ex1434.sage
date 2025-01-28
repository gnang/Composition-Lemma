# Loading the Hypermatrix Algebra Package.
load('Hypermatrix_Algebra_Package_code.sage')

# Initialization of the number of vertices
sz=Integer(5)

# Initialization of the function Tf
Tf=local_iterate_tuple([(Integer(0),Integer(0))]+\
[(i,i-Integer(1)) for i in rg(Integer(1),sz)], [sz-Integer(1)])

# Initialization of the local iterate of Tf that we call here Tg
Tg=local_iterate_tuple(Tf, tpl_pre_image_set(Tf,Tf[sz-1][1]))

# Displaying the functional tree on screen
# Tuple2DiGraphII(Tf,sz).plot().show(); Tuple2DiGraphII(Tg,sz).plot().show()

# Initialization of the list of vertex variables
X=var_list('x',sz)

# Initialization of the vertex Vandermonde factor
V=prod(X[j]-X[i] for i in rg(sz) for j in rg(sz) if i<j)

# Initialization of the factor which checks for overlaps between edge labels except for
# leaf edges, outgoing from sz-1 and its sibligns.
Ef1=prod((X[Tf[v][1]]-X[v]) + ((-Integer(1))^t)*(X[Tf[u][1]]-X[u]) for \
t in rg(2) for u in rg(1+Tf[sz-1][1]) for v in rg(1+Tf[sz-1][1]) if u<v)

# Initialization of the factor which checks for overlaps among leaf edges leaf edges,
# outgoing from sz-1 and its sibligns
Ef2=prod((X[Tf[v][1]]-X[v]) + ((-Integer(1))^t)*(X[Tf[u][1]]-X[u]) for \
t in rg(2) for u in rg(1+Tf[sz-1][1]) for v in tpl_pre_image_set(Tf, Tf[sz-1][1]))

# Initialization of the factor which checks for overlaps across leaf edges leaf edges,
# outgoing from sz-1 and its sibligns and other edges.
Ef3=prod((X[Tf[v][1]]-X[v]) + ((-Integer(1))^t)*(X[Tf[u][1]]-X[u]) for \
t in rg(2) for v in tpl_pre_image_set(Tf, Tf[sz-1][1]) for u in rg(v) if Tf[sz-1][1]<u)

# Initialization of the polynomial Pf
Pf=V*Ef1*Ef2*Ef3

# Initialization of the canonical representative i.e. the remainder
cPf=LagrangeInterpolationZn(Pf,X)

# Initialization of the factor which checks for overlaps between edge labels except for
# leaf edges, outgoing from sz-1 and its sibligns.
Eg1=prod((X[Tf[v][1]]-X[v]) + ((-Integer(1))^t)*(X[Tf[u][1]]-X[u]) for \
t in rg(2) for u in rg(1+Tf[sz-1][1]) for v in rg(1+Tf[sz-1][1]) if u<v)

# Initialization of the factor which checks for overlaps among leaf edges leaf edges,
# outgoing from sz-1 and its sibligns
Eg2=prod((X[Tg[v][1]]-X[v]) + ((-Integer(1))^t)*(X[Tf[u][1]]-X[u]) for \
t in rg(2) for u in rg(1+Tf[sz-1][1]) for v in tpl_pre_image_set(Tf, Tf[sz-1][1]))

# Initialization of the factor which checks for overlaps across leaf edges leaf edges,
# outgoing from sz-1 and its sibligns and other edges.
Eg3=prod((X[Tg[v][1]]-X[v]) + ((-Integer(1))^t)*(X[Tg[u][1]]-X[u]) for \
t in rg(2) for v in tpl_pre_image_set(Tf, Tf[sz-1][1]) for u in rg(v) if Tf[sz-1][1]<u)

# Initialization of the polynomial Pg
Pg=V*Eg1*Eg2*Eg3

# Initialization of the canonical representative i.e. the remainder
cPg=LagrangeInterpolationZn(Pg,X)


# ------ Telescoping Setup ------

# Initialization variables which helps keep track of the summand color scheme.
r,b=var('r,b')

# Initialization of the first factor
H1=prod((X[Tf[v][1]]-X[v]) + ((-Integer(1))^t)*(X[Tf[u][1]]-X[u]) \
for t in rg(2) for u in rg(1+Tf[sz-1][1]) for v in rg(1+Tf[sz-1][1]) if u<v)

# Initialization of the second factor
H2=prod((X[Tg[v][1]]-X[Tf[v][1]])*b + (X[Tf[v][1]]-X[v])*r + \
((-Integer(1))^t)*(X[Tf[u][1]]-X[u])*r for t in rg(2) for u in \
rg(1+Tf[sz-1][1]) for v in tpl_pre_image_set(Tf, Tf[sz-1][1])) 

# Initialization of the third factor A for which we neglect the color
H3a=prod((X[Tf[v][1]]-X[v])   - (X[Tf[u][1]]-X[u])   + \
0*(X[Tg[v][1]]-X[Tf[v][1]])   for v in \
tpl_pre_image_set(Tf, Tf[sz-1][1]) for u in rg(v) if Tf[sz-1][1]<u)
# Initialization of the third factor B which accounts for the colors
H3b=prod((X[Tf[v][1]]-X[v])*r + (X[Tf[u][1]]-X[u])*r + \
2*(X[Tg[v][1]]-X[Tf[v][1]])*b for v in \
tpl_pre_image_set(Tf, Tf[sz-1][1]) for u in rg(v) if Tf[sz-1][1]<u)

# Expanding the chromatic part
Hc=expand(H2*H3b)

# Initialization of the bichromatic part which removed which excludes
# the monochromatic red part.
HcFr=fast_reduce(Hc, [r^Hc.degree(r)], [0]).subs([b==Integer(1), r==Integer(1)])
Rfg=V*H1*H3a*HcFr


# Checking the telescoping
print('((V*H1*H2*H3a*H3b).subs([b==Integer(1), r==Integer(1)])-Pg).is_zero() evaluates to ',\
((V*H1*H2*H3a*H3b).subs([b==Integer(1), r==Integer(1)])-Pg).is_zero())

print('((V*H1*H2*H3a*H3b).subs([b==Integer(0), r==Integer(1)])-Pf).is_zero() evaluates to ',\
((V*H1*H2*H3a*H3b).subs([b==Integer(0), r==Integer(1)])-Pf).is_zero())

print('(Pf+Rfg-Pg).is_zero() evaluates to ',(Pf+Rfg-Pg).is_zero())

# Initialization of the canonical representative of Rfg
crRfg=LagrangeInterpolationZn(Rfg,X)

# Isolating summand which features the sub-monochromatic Blue part
Rfgb=V*H1*Hc.subs([b==Integer(1), r==Integer(0)])
# Obtaining the canonical representative or remainder
crRfgb=LagrangeInterpolationZn(Rfgb, X)

# Identifying the permutation which results in a graceful labeling of Gg
Lsig=[T for T in PermutationFunctionList(sz) if 0!=Pg.subs([X[i]==T[i][1] for i in rg(sz)])]

# Initialization of the symmetrized Rfg exploiting the common factor
scrRfg=sum(HcFr.subs([X[i]==T[i][1] for i in rg(sz)])*LagrangeBasis(T,X,sz) for T in Lsig)

# Isolating monochromatic Red part
Pfr=V*H1*H3a*(Hc.subs([b==Integer(0), r==Integer(1)]))
crPfr=LagrangeInterpolationZn(Pfr, X)

# Obtaining the automorphism group of the canonical representative of
# Rfg and its symmetrized version
AutGrp0=automorphism_group(crRfg,X); AutGrp1=automorphism_group(scrRfg,X)