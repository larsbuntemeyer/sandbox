#!/usr/bin/python

from employee import *

"This would create first object of Employee class"
emp1 = Employee("Zara", 2000)
"This would create second object of Employee class"
emp2 = Employee("Manni", 5000)

emp1.displayEmployee()
emp2.displayEmployee()
print "Total Employee %d" % Employee.empCount

emp1.age = 7  # Add an 'age' attribute.
emp1.age = 8  # Modify 'age' attribute.
#del emp1.age  # Delete 'age' attribute.

print hasattr(emp1, 'age')    # Returns true if 'age' attribute exists
print getattr(emp1, 'age')    # Returns value of 'age' attribute
print setattr(emp1, 'age', 8) # Set attribute 'age' at 8
print delattr(emp1, 'age')    # Delete attribute 'age'

print "Employee.__doc__:", Employee.__doc__
print "Employee.__name__:", Employee.__name__
print "Employee.__module__:", Employee.__module__
print "Employee.__bases__:", Employee.__bases__
print "Employee.__dict__:", Employee.__dict__
