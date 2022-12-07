function dtstring = DT4filename

dtdt=datetime;
         dty =yyyymmdd(dtdt); dty = num2str(dty);dty= dty(3:end);
        [h,m,s] = hms(dtdt);
                hh = sprintf('%02d',h);
                mm = sprintf('%02d',m);
                ss = sprintf('%02d',round(s));
                dtstring = strcat('_DT',dty,'-',hh,mm,ss);

end