# Name:- VIKRANT CHAUHAN
# Roll No. :- 231110407
# Batch:- Y23
# Programe:- MS(R) CYBERSECURITY
import math
import numpy as np

n=512
q=12289
gamma=10968
omega_HH=(gamma**2)%q

gamma_array=[0]*n
gamma_inverse_array=[0]*n


def make_gamma_array():
    for i in range(0,n):
        temp=pow(gamma,i,q)
        gamma_array[i]=temp



def n_inverse(number):
    for i in range(0,q):
        if(((i*number)%q)==1):
            return i
        

def gamma_inverse_calculate(g):
    for i in range(2,q):
        if((i*g)%q==1):
            return i
        

def reverse_num_binary(number):
    bits=int(math.log2(n))
    reversed_number = 0
    for i in range(bits):
        reversed_number <<= 1
        reversed_number |= (number & 1)
        number >>= 1
    return reversed_number


def make_gamma_inverse_array():
    gamma_inv=gamma_inverse_calculate(gamma)
    for i in range(0,n):
        temp=pow(gamma_inv,i,q)
        t=reverse_num_binary(i)
        gamma_inverse_array[t]=temp


def Bitreverse(original_list,bitSize):
    modified_list=[]
    for i in range(0,n):
        binary_num= format(int(i), 'b')
        reversed_binary_num=binary_num
        if(len(binary_num)<bitSize):
            temp_str="0"*int(bitSize-len(binary_num))
            reversed_binary_num=temp_str+binary_num
        reversed_binary_num = reversed_binary_num[::-1]
        modified_list.append(original_list[int(reversed_binary_num,2)])
    return modified_list


def NTT(P_NN):
    P_temp_Modified=P_NN
    t = n
    m=1
    gamma_original_ntt=Bitreverse(gamma_array,math.log2(n))
    while(m<n):
        t = int(t/2)
        for i in range(0,m):
            j1 =2*i*t
            j2 = j1 + t - 1
            S = gamma_original_ntt[m + i]
            for j in range(j1,j2+1):
                u1 = P_temp_Modified[j]
                t1 =P_temp_Modified[j+t]*S 
                P_temp_Modified[j]=(u1+t1)%q 
                P_temp_Modified[j+t]=(u1-t1)%q
        m=m*2
    return P_temp_Modified


def Inverse_NTT(P_Modified):
    P_temp_Modified2=P_Modified
    t = 1
    m=n
    while(m>1):
        j1 = 0
        h = int(m/2)
        for i in range(0,h):
            j2 = j1 + t- 1
            G = gamma_inverse_array[h + i]
            j=j1
            
            while(j<=j2):
                u1 = P_temp_Modified2[j]
                t1 =P_temp_Modified2[j+t]
                P_temp_Modified2[j]=(u1+t1)%q
                P_temp_Modified2[j + t] = ((u1-t1)*G)%q   
                j+=1 
            j1 = j1 + 2*t 
        t = 2*t 
        m=int(m/2)
    for i in range(0,n):
        P_temp_Modified2[i]=(P_temp_Modified2[i]*n_inverse(n))%q
    return P_temp_Modified2


def PointWise_Multiplication(NTT_Poly1,NTT_Poly2):
    result=[]
    for i in range(0,len(NTT_Poly1)):
        temp=(NTT_Poly1[i]*NTT_Poly2[i])%q
        result.append(temp)
    return result


P_N1=np.random.randint(0, q, n)
print("PN1 is ",P_N1)
P_N2=np.random.randint(0, q, n)
print("PN2 is ",P_N2)

t1=[0]*n
t2=[0]*n

for i in range(0,n):
    t1[i]=P_N1[i]
    t2[i]=P_N2[i]

make_gamma_array()
print("Gamma array ",gamma_array)

P_N1_Modified1=NTT(t1)
print("Result Of NTT PN1 ",P_N1_Modified1)

P_N2_Modified1=NTT(t2)
print("Result Of NTT PN2 ",P_N2_Modified1)


Product_Result=PointWise_Multiplication(P_N1_Modified1,P_N2_Modified1)
print("Result Of Point wise Multiplication",Product_Result)

make_gamma_inverse_array()
print("Gamma inverse array ",gamma_inverse_array)

P_Final=Inverse_NTT(Product_Result)
print("Result Of INTT Final Answer ",P_Final)

print("Verification :)")
print("PN1 is ",P_N2)
print("PN2 is ",P_N2)

f = np.zeros(n + 1)
f[0] = 1
f[n] = 1
ans=np.remainder(np.polydiv(np.polymul(P_N1[::-1], P_N2[::-1]),f)[1],q).astype(int)[::-1]

print("Polynomial Multiplication and modulo ",ans)
print(np.array_equal(ans, P_Final))

if(np.array_equal(ans, P_Final)):
    print("Your Answer Is Correct")
else:
    print("Wrong")

