function dqi_dvecY = dqi_dvecY(VecX, VecY, MOD)
    [~, dqi_dvecY] = fqei_dfqeiXYU(VecX, VecY, [], MOD, 'q', 'i', 'Y');
end
