function remove_blankline,header

; INPUT:
;       Header = String array containing FITS header
;
; COMMENT:
;     Old Version: The routine assumes blank lines do exist in the header
;Version May 2018: The routine DOESN'T assume blank lines do exist - Castellanos Duran

blank_str=string(bytarr(80)+32B)
zblank=where(header eq blank_str,nblank)

nheader=n_elements(header)

if nblank gt 0 then begin
if(zblank[nblank-1] eq nheader-1) then begin
   header=header[0:nheader-2]
endif else begin
   header=[header[0:zblank[nblank-1]-1],header[zblank[nblank-1]+1:*]]
endelse
endif

return,header
end
