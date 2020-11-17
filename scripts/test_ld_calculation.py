from __future__ import division

p = 0.2
q = 0.5

u_1 = 0.1
u_2 = 0.3
v_1 = 0.5
v_2 = 0.2

k = u_1+ u_2 + v_1 + v_2

partial_f_p = q- (q**2) - 2*p*q + 2*p*(q**2)
partial_f_q = p- (p**2) - 2*p*q + 2*q*(p**2)

X = p*q*(1-p) * (1-q)

N = 1000000

def final_form():

    return -1*(2*N*k)*X + N*v_1*q*(1-q) + N*v_2*p*(1-p) + N*(u_1-v_1)*p*q*(1-q) + N*(u_2-v_2)*p*q*(1-p)


def initial_form():

    return N*(v_1 - (u_1+v_1)*p ) *partial_f_p +  N*(v_2 - (u_2+v_2)*q ) *partial_f_q


print(final_form() == initial_form())

print(final_form()
)
print(initial_form())
