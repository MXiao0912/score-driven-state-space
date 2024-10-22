c = [log(0.7);log(0.5);log(0.1);zeros(2,1)];
A = InvMaxZero(1,zeros(5,1)); 
B = InvMaxZero(0.1, zeros(5,1));

Init_tvp = [InvUnoMenoUno(psi);InvUnoMenoUno(phi);0;0;0;kap_hes;c;A;B];
Init_base = [InvUnoMenoUno(psi),InvUnoMenoUno(phi),0,0,0,exp(c(1:3))'];


base_kf_original(Init_base,y)

kf_tvp(Init_tvp,y, c)
