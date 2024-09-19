import numpy

a = numpy.arange(12).reshape(3,4)
print "a =", a

b = numpy.array([0.1, 0.2, 0.3, 0.4])
print "b =", b

ab = numpy.dot(a,b)
print "ab = dot(a,b) =", ab

c = numpy.array([0.1, 0.2, 0.3])
print "c =", c

print "ab + c =", ab + c
print "dot(ab,c) =", numpy.dot(ab,c)

d = a[0:3,0:3]
print "d = a[0:3,0:3] =", d

f = a[:,3]
print "f = a[:,3] =", f
