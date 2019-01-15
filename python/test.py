from __future__ import print_function
import numpy as np
a=np.array([2.1,3.2])
b=np.array([1.2,0.1])
import volesti
p=volesti.Point(a)
p2=volesti.Point(b)
u=np.array(p+p2)
a=np.array([[1,0],[0,1],[-1,0], [0,-1]])
b=np.array([1,1,1,1])
#b.shape = (2, 1);
p=volesti.Polytope(a, b)
print('Intersection tests')

p.create_point_representation(np.zeros(2))
source=np.array([0,0])
direction=np.array([1,0.5])
direction=direction/np.linalg.norm(direction);
print('line intersect', p.line_intersect(source, direction))
print('intersection with 0 facet', p.intersect(source, direction, 0))
print('boundary oracle', p.boundary(source, direction, 0.01, 100));

print('Begin membership tests')
dim=10
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

print('Begin boundary tests')
from timeit import default_timer as timer
failures = 0
boundary_time = 0.0
li_time = 0.0
avg_dists=0.0
for i in range(1000):
    direction=np.random.multivariate_normal(mean, cov, 1)
    direction=direction/np.linalg.norm(direction)
    source = np.zeros(dim)
    tic = timer()
    p1=p.boundary(source, direction, 0.001, 100) 
    p2=p.boundary(source, -direction, 0.001, 100) 
    toc = timer()
    boundary_time += toc - tic
    tic = timer()
    points = p.line_intersect(source, direction)
    toc = timer()
    li_time += toc -tic

    avg_dists += np.linalg.norm(p1-points[0,:])
    avg_dists += np.linalg.norm(p2-points[1,:])

    if np.linalg.norm(p1-points[0,:])>0.01:
        failures += 1
    if np.linalg.norm(p2-points[1,:])>0.01:
        failures += 1

print('Boundary time: {}.\nLI time: {}\nFailures:{}\nAvg dists:{}'.format(boundary_time/1000, li_time/1000, failures, avg_dists/2000))
