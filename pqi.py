import math
# Full IJS

def poly_coeff_from_roots(roots):
    coeffs = []
    
    # Go through all the roots
    n = len(roots)
    #print(n)
    
    # Generate the result Polynomial, which initially is set to -roots[0] + x :: index of the list can be seen as power of x
    coeffs.insert(0, (-1)*roots[0])
    coeffs.insert(1,1)
    
    # Continue multiplying -roots[i] + x
    i = 1
    while i < n:
        # right shift, equivalent to multiplying x, the index increased 1 as the result
        coeffs.insert(0,0)        
        j = 0
        while j < (len(coeffs) - 1):
            # mu
            coeffs[j] = coeffs[j] + coeffs[j+1] * (-1) * roots[i]
            j = j + 1

        i = i + 1
    return coeffs[::-1]  ## reverse list from highest power coefficient to lowest

#print(poly_coeff_from_roots([1,2,3,4,5]))
##################################################################################################################
def get_coeff_list(c, w, t, v):
    """
    c - sent coeffs
    w - secondary list at second device

    s, n,
    s - size of w ##### ## b
    n - window size for root selection (not used -- same as s) #### number of roots to generate s-(t/2)

    v - number of varations of roots 
    ### IMPORTANT: CAN ONLY CORRECT ONE ERROR or t/2 errors
    """
    ##################
    s = len(w)
    n = s 
    #print(c, w, s, n, t, v)
    ##################

    result = [] # coefficients found
    corr_seq = []
    #print(c)
    # # For range finding on second device (not used if window is fixed or full)
    #for i in range(0, s-n+1, 1):
        #r_w = w[i: i+n] # values of w in the window - to find root of... (should be full - in algo.)
    ### Indent if using range -- till before return
    r_w = w.copy()
    #print(r_w)
    for err_idx in range(n):
        removed_r = r_w[err_idx]
        k = int(math.pow((2 * v + 1), (n-1)))
        #print(k)
        while k > 0:
            temp_list = r_w.copy()
            temp_list.remove(removed_r)
            temp_variation = k - 1

            for l in range(n-1):
                temp_list[l] = temp_list[l] + (temp_variation % (2 * v + 1)) - v
                temp_variation = temp_variation // (2 * v + 1)
            #print(temp_list)
            residue_coeff = poly_coeff_from_roots(temp_list)
            #print(residue_coeff)
            addd_rt = (-1) * (c[1] - residue_coeff[1])
            temp_list.insert(err_idx,addd_rt)
            residue_coeff = poly_coeff_from_roots(temp_list)
            #print(residue_coeff)
            k = k - 1
            if(residue_coeff not in result):
                #print(residue_coeff)
                coeff_check = True
                
                for i_t in range(len(c)-1): #+1 -1 # # t
                    if(c[i_t+1] != residue_coeff[i_t+1]): # or not coeff_check
                        coeff_check = False
                        break
                if(coeff_check):
                    corr_seq = temp_list
                    result.append(residue_coeff)
                # if(c[1] == residue_coeff[1] and c[2] == residue_coeff[2]):#  and c[3] == residue_coeff[3]):
                #     result.append(residue_coeff)
        ###
    print("Corrected Sequence at Bob:", corr_seq)
    if(len(result)>1):
        print("Multiple Sequences")
    return result


################################################
################################################
######### SPLINE

import numpy as np

def spline_projection(seq):
    projected_seq = []
    for i in range(1, len(seq)-1, 1):
        A = np.array(
            [[(i-1)**3,     (i-1)**2,   (i-1)**1,   1,  0,       0,      0,      0],
            [0,             0,          0,          0,  i**3,    i**2,   i**1,   1],
            [i**3,          i**2,       i**1,       1,  0,       0,      0,      0],
            [(i+1)**3,      (i+1)**2,   (i+1)**1,   1,  0,       0,      0,      0],
            [3*i**2,        2*i,        1,          0,  -3*i**2, -2*i,   -1,     0],
            [6*i,           2,          0,          0,  -6*i,    -2,     0,      0],
            [6*(i-1),       2,          0,          0,  0,       0,      0,      0],
            [0,             0,          0,          0,  6*(i+1), 2,      0,      0]]
        )
        b = np.array([seq[i-1], seq[i], seq[i], seq[i+1], 0, 0, 0, 0])
        coeff = np.dot(np.linalg.inv(A), b)
        
        #print(np.linalg.det(A))
        #print(np.dot(np.linalg.inv(A), b))
        #print(i)

        k = i+2
        # projected cubic eq.
        projected_seq.append( round(coeff[0]*k**3 + coeff[1]*k**2 + coeff[2]*k**1 + coeff[3]) )
        #projected_seq.append( round(coeff[4]*k**3 + coeff[5]*k**2 + coeff[6]*k**1 + coeff[7]) )
    return projected_seq

# # TESTS
# print(spline_projection(alice))
# print()
# print(spline_projection(bob))

################################################
################################################
#Quantization - uniform

# in case the range is too low
global previous_lvls
previous_lvls = []
def quantization_m(data, m):
    global previous_lvls
    
    data = np.array(data)
    quant_seq = []

    variance = (1/len(data)) * np.sum(data**2)
    norm_x = data/np.sqrt(variance)
    
    x_min = min(norm_x)
    x_max = max(norm_x)

    delta = (x_max-x_min)/m

    for x in norm_x:
        # lvls = np.arange(x_min, x_max+1, delta)

        try:
            lvls = np.arange(x_min, x_max+1, delta)
            #print(lvls)
        except:
            lvls = previous_lvls

        previous_lvls = lvls

        for i in range(1,len(lvls)):
            if(x>=lvls[i-1] and x<=lvls[i]):
                #print(i)
                quant_seq.append(i)
                break
    return quant_seq
