function dqe_dvecX = dqe_dvecX(MOD, VecX, VecY, ~)
    [~, dqe_dvecX] = MOD.fqei_dfqeiXYU([VecX; VecY], 'q', 'e', 'X');
end
