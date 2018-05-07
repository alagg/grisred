function remove_comment,header

; INPUT:
;       Header = String array containing FITS header
;


comment_str='COMMENT   '+string(bytarr(70)+32B)
zcomment=where(header eq comment_str,ncomment)

nheader=n_elements(header)

if ncomment gt 0 then begin
if(zcomment[ncomment-1] eq nheader-1) then begin
   header=header[0:nheader-2]
endif else begin
   header=[header[0:zcomment[ncomment-1]-1],header[zcomment[ncomment-1]+1:*]]
endelse
endif

return,header
end
