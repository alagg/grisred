;;
;; Routine to add new keyword such as the GIT revision keywords to the header,
;; eliminate unnecessary whitespaces, and the add at the end of the header 
;; "END" (standard)
;;
;; V1. May 2018 - Sebastian Castellanos Duran

pro patch_hdr_v6,map_out,keyw,keyv,keyc

data = readfits(map_out+'cc',hdrstr)

print,'Writing keywords on header...'
for j=0,n_elements(keyw) -1 do begin

  print,keyw[j]+' -> '+keyv[j]+' -> '+keyc[j]
  tmp = keyv[j]
  if is_number2(tmp) then tmp = double(keyv[j])
  sxaddpar,hdrstr,keyw[j],tmp,keyc[j],before='LC1-1'

endfor

;;remove blank spaces
blank_str=string(bytarr(80)+32B)
zblank=where(hdrstr eq blank_str,nblank)
if nblank gt 0 then for j=0,nblank-1 do hdrstr=remove_blankline(hdrstr)

WRITEFITS, map_out+'cc', data, hdrstr
end
