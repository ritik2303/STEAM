function dqi_dvecX = dqi_dvecX(VecX, VecY, MOD)
    [~, dqi_dvecX] = fqei_dfqeiXYU(VecX, VecY, [], MOD, 'q', 'i', 'X');
end
