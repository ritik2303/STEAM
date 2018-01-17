function dfe_dvecX = dfe_dvecX(MOD, VecX, VecY, VecU, ~)
    [~, dfe_dvecX] = MOD.fqei_dfqeiXYU([VecX; VecY], 'f', 'e', 'X');
end
