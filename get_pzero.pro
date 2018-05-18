function get_pzero,date
  
  year=fix(strmid(date,0,4))
  month=fix(strmid(date,5,2))
  day=fix(strmid(date,8,2))
  
  
  case 1 of:
    year lt 2016: pzero=-1.2
    year eq 2016: pzero=53.6
    year eq 2017 and month le 4: pzero=53.6
    year eq 2017 and month gt 4 and month le 7: pzero=41.8 ; May 1st, 2017 
    year eq 2017 and month gt 7: pzero=53.6                ; Aug 31st, 2017 
    year eq 2018: pzero=-1.8                               ; Apr 13st, 2018 
    else: begin
      pzero=-1.8
      message,/cont,'Unknown time for pzero computation ('+date+')'
      message,/cont,'set pzero='+string(pzero)
    endelse
  endcase
  
  print,'PZERO = '+string(pzero,format='(f6.2)')+' ('+date+')'
  
  stop
  return,pzero
end
