# Macro for Park transformation*/
AB2M=lambda A, B, COS, SIN : ( (A)*COS  + (B)*SIN )
AB2T=lambda A, B, COS, SIN : ( (A)*-SIN + (B)*COS ) 
MT2A=lambda M, T, COS, SIN : ( (M)*COS - (T)*SIN )
MT2B=lambda M, T, COS, SIN : ( (M)*SIN + (T)*COS )

isNumber = lambda x : isinstance(x,(np.int,np.int64,np.float,np.float64))

# Macro for Power-invariant inverse Clarke transformation
AB2U=lambda A, B : ( 0.816496580927726 * ( A ) )
AB2V=lambda A, B : ( 0.816496580927726 * ( A*-0.5 + B*0.8660254037844387 ) )
AB2W=lambda A, B : ( 0.816496580927726 * ( A*-0.5 + B*-0.8660254037844385 ) )