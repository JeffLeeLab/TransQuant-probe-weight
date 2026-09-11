function ims2 = enhance_image(ims,method)
method
if strcmp(method,'histeq')
    ims2=histeq(uint16(ims));
elseif strcmp(method,'adapthisteq')
    ims2=adapthisteq(uint16(ims));
else strcmp(method,'imadjust')
    ims2=imadjust(uint16(ims));
end
    