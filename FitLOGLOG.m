function FitLOGLOG(x,y,ordre)
    P = polyfit(log(x),log(y),ordre); 
    z  = polyval(P, log(x));
    loglog(x,exp(z),'--','Linewidth',1);
    legendStrings = string(P(1));
    leg = legend(legendStrings,'Location','northeast','NumColumns',2);
    title(leg, 'linear fit: slope')
end