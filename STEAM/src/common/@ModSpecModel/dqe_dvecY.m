function dqe_dvecY = dqe_dvecY(MOD, VecX, VecY, ~)
    [~, dqe_dvecY] = MOD.fqei_dfqeiXYU([VecX; VecY], 'q', 'e', 'Y');
end
