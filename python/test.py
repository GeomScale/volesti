import numpy as np
a=np.array([2.1,3.2])
b=np.array([1.2,0.1])
import volesti
p=volesti.Point(a)
p2=volesti.Point(b)
u=np.array(p+p2)
print(np.array(p2))
print(u)

a=np.array([[2,0],[0,2],[-1,0], [0,-1]])
b=np.array([2,2,1,1])
#b.shape = (2, 1);
p=volesti.Polytope(a, b)
print('A=', a)
print('b=', b)

print('Begin membership tests')
dim=100
n=100000
mean=np.zeros(dim)
cov=np.eye(dim)
normals=np.random.multivariate_normal(mean, cov, n)
b=np.ones(n)
p=volesti.Polytope(normals, b)
p.create_point_representation(mean)

test_point = np.zeros(dim)
print(p.contains_point(test_point), p.is_in(test_point))
test_point[0]=1.0
print(p.contains_point(test_point), p.is_in(test_point))
