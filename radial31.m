function retval =  radial31(mat,a)

retval = (1-a*mat);
retval(retval<0)=0;
retval = retval.^4.*(4*mat+1);

end