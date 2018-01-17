function fe = fe(MOD, VecX, VecY, VecU, ~)
    fe = MOD.fqei_dfqeiXYU([VecX; VecY], 'f', 'e');
end
