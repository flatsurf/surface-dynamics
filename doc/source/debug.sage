print("L0: "); from surface_dynamics import iet
print("L1: "); perm = iet.Permutation('a b c d', 'd c b a')
print("L2: "); x = polygen(QQ)
print("L3: "); K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
print("L4: "); length2 = [1, sqrt2, sqrt2**2, sqrt2**3]
print("L5: "); t2 = iet.IntervalExchangeTransformation(perm, length2)
print("L6: "); print(t2)
print("L7: "); from surface_dynamics.interval_exchanges.conversion import iet_to_pyintervalxt, iet_from_pyintervalxt
print("L8: "); u2 = iet_to_pyintervalxt(t2)
print("L9: "); print(u2)
print("L10: "); v2 = iet_from_pyintervalxt(u2)
print("L11: "); print(v2)
print("L12: "); print(u2.boshernitzanNoPeriodicTrajectory())
