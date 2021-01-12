import numpy as np
import matplotlib.pyplot as plt



def GDA(func_w, func_a, eta_c, w, alpha, maxiter,c):
    w_list = []
    w_list.append(w)
    alpha_list = []
    alpha_list.append(alpha)
    beta = alpha
    for t in range(maxiter):
        fg_w = func_w(w,alpha)
        fg_alpha = func_a(w,alpha)
        
        if eta_c =='1/t':
            eta = 1/(t+1)
        elif eta_c =='1/t^2':
            eta = 1/(t+1)**2
        elif eta_c == 'sq':
            eta = 1/(t+1)**(1/2)
        
        w = w-eta*fg_w
        w_list.append(w)
        
        N = alpha + eta*fg_alpha
        alpha = Proj(N,c)
        alpha_list.append(alpha)
        
    return w_list, alpha_list



def Proj(N,c):
    result = []
    for i in N:
        result.append(min(10, max(0,i)))
    return np.array(result)



def func_w(w, alpha):
    result = np.zeros(len(w))
    for i in range(len(alpha)):
        result += y[i]*X[i]*alpha[i]
    return w - result



def func_alpha(w, alpha):
    result = []
    for i in range(len(alpha)):
        result.append(1- np.dot(y[i]*X[i],w))  
    return np.array(result)



def draw(w1,a1, s):
    w_0 = [item[0] for item in w1]
    w_1 = [item[1] for item in w1]
    w_2 = [item[2] for item in w1]
    a_0 = [item[0] for item in a1]
    a_1 = [item[1] for item in a1]
    a_2 = [item[2] for item in a1]
    a_3 = [item[3] for item in a1]
    plt.plot(range(1001),w_0,label = 'w[0]')
    plt.plot(range(1001),w_1,label = 'w[1]')
    plt.plot(range(1001),w_2,label = 'w[2]')
    plt.legend()
    if(s == 1):
        plt.title('w vs. iteration 1/t')
    elif(s == 2):
        plt.title('w vs. iteration 1/t^2')
    else:
        plt.title('w vs. iteration 1/sq(t)')


    plt.show()
    
    plt.plot(range(1001),a_0,label = 'a[0]')
    plt.plot(range(1001),a_1,label = 'a[1]')
    plt.plot(range(1001),a_2,label = 'a[2]')
    plt.plot(range(1001),a_3,label = 'a[3]')
    plt.legend()
    if(s == 1):
        plt.title('a vs. iteration 1/t')
    elif(s == 2):
        plt.title('a vs. iteration 1/t^2')
    else:
        plt.title('a vs. iteration 1/sq(t)')
    plt.show()



def Avg_GDA(func_w, func_a, eta_c, w, alpha, maxiter,c):
    w_list = []
    w_list.append(w)
    alpha_list = []
    alpha_list.append(alpha)
    beta = alpha
    sum_eta = 0
    sum_w = 0
    sum_a = 0
    for t in range(maxiter):
        fg_w = func_w(w,alpha)
        fg_alpha = func_a(w,alpha)
        
        if eta_c =='1/t':
            eta = 1/(t+1)
        elif eta_c =='1/t^2':
            eta = 1/(t+1)**2
        elif eta_c == 'sq':
            eta = 1/(t+1)**(1/2)
        
        sum_eta += eta
        
        w = w-eta*fg_w
        
        
        N = alpha + eta*fg_alpha
        alpha = Proj(N,c)
        sum_w += eta*w
        sum_a += eta*alpha
        
        w_avg = sum_w/sum_eta
        a_avg = sum_a/sum_eta
        w_list.append(w_avg)
        alpha_list.append(a_avg)
        
        
        
    return w_list, alpha_list




def draw2(w1,a1, s):
    w_0 = [item[0] for item in w1]
    w_1 = [item[1] for item in w1]
    w_2 = [item[2] for item in w1]
    a_0 = [item[0] for item in a1]
    a_1 = [item[1] for item in a1]
    a_2 = [item[2] for item in a1]
    a_3 = [item[3] for item in a1]
    plt.plot(range(1001),w_0,label = 'w[0]')
    plt.plot(range(1001),w_1,label = 'w[1]')
    plt.plot(range(1001),w_2,label = 'w[2]')
    plt.legend()
    if(s == 1):
        plt.title('#2.4 w vs. iteration 1/t')
    elif(s == 2):
        plt.title('#2.4 w vs. iteration 1/t^2')
    else:
        plt.title('#2,4 w vs. iteration 1/sq(t)')


    plt.show()
    
    plt.plot(range(1001),a_0,label = 'a[0]')
    plt.plot(range(1001),a_1,label = 'a[1]')
    plt.plot(range(1001),a_2,label = 'a[2]')
    plt.plot(range(1001),a_3,label = 'a[3]')
    plt.legend()
    if(s == 1):
        plt.title('#2.4 a vs. iteration 1/t')
    elif(s == 2):
        plt.title('#2.4 a vs. iteration 1/t^2')
    else:
        plt.title('#2.4 a vs. iteration 1/sq(t)')
    plt.show()



def Nest_GDA(func_w, func_a, eta_c, w, a, maxiter,c,beta):
    w_list = []
    w_list.append(w)
    alpha_list = []
    alpha_list.append(a)
    prev_w = w
    prev_a = a
    for t in range(maxiter):
        
        
        if eta_c =='1/t':
            eta = 1/(t+1)
        
        '''
        if t == 0:
            continue
        '''
        
        w_tuta = w + beta*(w-prev_w)
        a_tuta = a + beta*(a-prev_a)
        prev_w = w
        prev_a = a
        
        fg_w = func_w(w_tuta,a_tuta)
        fg_alpha = func_a(w_tuta,a_tuta)
        w = w_tuta-eta*fg_w
        
        
        N = a_tuta + eta*fg_alpha
        a = Proj(N,c)
       
        
        
        w_list.append(w)
        alpha_list.append(a)
        
        
        
    return w_list, alpha_list



def draw3(w1,a1, s):
    w_0 = [item[0] for item in w1]
    w_1 = [item[1] for item in w1]
    w_2 = [item[2] for item in w1]
    a_0 = [item[0] for item in a1]
    a_1 = [item[1] for item in a1]
    a_2 = [item[2] for item in a1]
    a_3 = [item[3] for item in a1]
    plt.plot(range(1001),w_0,label = 'w[0]')
    plt.plot(range(1001),w_1,label = 'w[1]')
    plt.plot(range(1001),w_2,label = 'w[2]')
    plt.legend()
    if(s == 1):
        plt.title('#2.5 w vs. beta 0.1')
    elif(s == 2):
        plt.title('#2.5 w vs. beta -0.1')
    else:
        plt.title('#2.5 w vs. iteration 1/sq(t)')


    plt.show()
    
    plt.plot(range(1001),a_0,label = 'a[0]')
    plt.plot(range(1001),a_1,label = 'a[1]')
    plt.plot(range(1001),a_2,label = 'a[2]')
    plt.plot(range(1001),a_3,label = 'a[3]')
    plt.legend()
    if(s == 1):
        plt.title('#2.5 a vs. beta 0.1')
    elif(s == 2):
        plt.title('#2.5 a vs. beta -0.1')
    else:
        plt.title('#2.5 a vs. iteration 1/sq(t)')
    plt.show()



def Extra_GDA(func_w, func_a, eta_c, w, a, maxiter,c,beta):
    w_list = []
    w_list.append(w)
    alpha_list = []
    alpha_list.append(a)

    for t in range(maxiter):
        
        
        if eta_c =='1/t':
            eta = 1/37
        else:
            eta = 1
            p = 0
            fg_w = func_w(w,a)
            fg_alpha = func_a(w,a)
            
            while(2*eta > p):
                
                
                w_tuta = w - eta*fg_w
                N = a + eta*fg_alpha
                a_tuta = Proj(N,c)
                
                eta  = eta/2
                
                fg_w2 = func_w(w_tuta,a_tuta)
                fg_alpha2 = func_a(w_tuta,a_tuta)
        

                upper = np.linalg.norm(w - w_tuta)**2 + np.linalg.norm(a - a_tuta)**2
                lower = np.linalg.norm(fg_w - fg_w2)**2 + np.linalg.norm(fg_alpha - fg_alpha2)**2
                
                upper2 = upper**(1/2)
                lower2 = lower**(1/2)
                
                
                
                p = upper2/lower2/2
                
        
        
        fg_w = func_w(w,a)
        fg_alpha = func_a(w,a)
        
        w_tuta = w -eta*fg_w
        
        
        N = a + eta*fg_alpha
        a_tuta = Proj(N,c)
        
        fg_w2 = func_w(w_tuta,a_tuta)
        fg_alpha2 = func_a(w_tuta,a_tuta)
        
        w = w - eta*fg_w2
        
        N2 = a + eta*fg_alpha2
        a = Proj(N,c)
        
        
        w_list.append(w)
        alpha_list.append(a)
        
        
        
    return w_list, alpha_list



def draw4(w1,a1, s):
    w_0 = [item[0] for item in w1]
    w_1 = [item[1] for item in w1]
    w_2 = [item[2] for item in w1]
    a_0 = [item[0] for item in a1]
    a_1 = [item[1] for item in a1]
    a_2 = [item[2] for item in a1]
    a_3 = [item[3] for item in a1]
    plt.plot(range(1001),w_0,label = 'w[0]')
    plt.plot(range(1001),w_1,label = 'w[1]')
    plt.plot(range(1001),w_2,label = 'w[2]')
    plt.legend()
    if(s == 1):
        plt.title('#2.6 w vs. eta = 1/37')
    elif(s == 2):
        plt.title('#2.6 w vs. line search')
    else:
        plt.title('#2.5 w vs. iteration 1/sq(t)')


    plt.show()
    
    plt.plot(range(1001),a_0,label = 'a[0]')
    plt.plot(range(1001),a_1,label = 'a[1]')
    plt.plot(range(1001),a_2,label = 'a[2]')
    plt.plot(range(1001),a_3,label = 'a[3]')
    plt.legend()
    if(s == 1):
        plt.title('#2.6 a vs. eta = 1/37')
    elif(s == 2):
        plt.title('#2.6 a vs. line search')
    else:
        plt.title('#2.6 a vs. iteration 1/sq(t)')
    plt.show()


if __name__ == "__main__":
    X = np.array([[1,1],[1,-1],[-1,1],[-1,-1]])
    X = np.append(X,np.ones((4,1)),axis=1)
    y = np.array([1,1,1,-1])
    init_w = np.zeros(3)
    init_alpha = np.zeros(4)

    #2.3
    w1, a1 = GDA(func_w, func_alpha, '1/t', init_w, init_alpha, 1000,10)
    draw(w1,a1,1)

    w1, a1 = GDA(func_w, func_alpha, '1/t^2', init_w, init_alpha, 1000,10)
    draw(w1,a1,2)

    w1, a1 = GDA(func_w, func_alpha, 'sq', init_w, init_alpha, 1000,10)
    draw(w1,a1,3)

    #2.4

    w1, a1 = Avg_GDA(func_w, func_alpha, '1/t', init_w, init_alpha, 1000,10)
    draw2(w1,a1,1)

    w1, a1 = Avg_GDA(func_w, func_alpha, '1/t^2', init_w, init_alpha, 1000,10)
    draw2(w1,a1,2)

    w1, a1 = Avg_GDA(func_w, func_alpha, 'sq', init_w, init_alpha, 1000,10)
    draw(w1,a1,3)

    #2.5
    w1, a1 = Nest_GDA(func_w, func_alpha, '1/t', init_w, init_alpha, 1000,10,0.1)
    draw3(w1,a1,1)

    w1, a1 = Nest_GDA(func_w, func_alpha, '1/t', init_w, init_alpha, 1000,10,-0.1)
    draw3(w1,a1,2)

    #2.6

    w1, a1 = Extra_GDA(func_w, func_alpha, '1/t', init_w, init_alpha, 1000,10,-0.1)
    draw4(w1,a1,1)

    w1, a1 = Extra_GDA(func_w, func_alpha, 'line search', init_w, init_alpha, 1000,10,-0.1)
    draw4(w1,a1,2)




