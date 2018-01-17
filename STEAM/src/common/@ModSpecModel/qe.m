function qe = qe(MOD, VecX, VecY, ~)
    qe = MOD.fqei_dfqeiXYU([VecX; VecY], 'q', 'e');
end
