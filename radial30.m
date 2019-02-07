function retval =  radial30(mat,a)

retval = (1-a*mat);
retval(retval<0)=0;
retval = retval.^2;

end