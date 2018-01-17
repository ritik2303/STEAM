function dfi_dvecX = dfi_dvecX(VecX, VecY, VecU, MOD)
    [~, dfi_dvecX] = fqei_dfqeiXYU(VecX, VecY, [], MOD, 'f', 'i', 'X');
end
