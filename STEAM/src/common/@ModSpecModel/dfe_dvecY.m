function dfe_dvecY = dfe_dvecY(MOD, VecX, VecY, VecU, ~)
    [~, dfe_dvecY] = MOD.fqei_dfqeiXYU([VecX; VecY], 'f', 'e', 'Y');
end
