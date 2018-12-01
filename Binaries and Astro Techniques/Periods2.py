def orbitalPhase(P,To,t):

    P=float(P)
    To=float(To)
    t=float(t)
    ##This converts all values to double floating precision values  
    phi =((t-To))/P    
    if phi < 0:
        phi += 1
    while phi > 1:        
        
        n=2
        To = To + (P*n)
        phi = abs((t-To)/P)
        n +=1
        if phi%1==0.0:
            break
        
        
        
        
    return phi


