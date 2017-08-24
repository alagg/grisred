function remove_blankline,header

; INPUT:
;       Header = String array containing FITS header
;
; COMMENT:
;   The routine assumes blank lines do exist in the header

blank_str=string(bytarr(80)+32B)
zblank=where(header eq blank_str,nblank)

nheader=n_elements(header)

if(zblank[nblank-1] eq nheader-1) then begin
   header=header[0:nheader-2]
endif else begin
   header=[header[0:zblank[nblank-1]-1],header[zblank[nblank-1]+1:*]]
endelse

return,header
end
