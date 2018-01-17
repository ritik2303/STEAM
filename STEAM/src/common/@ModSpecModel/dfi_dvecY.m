function dfi_dvecY = dfi_dvecY(VecX, VecY, VecU, MOD)
    [~, dfi_dvecY] = fqei_dfqeiXYU(VecX, VecY, [], MOD, 'f', 'i', 'Y');
end
