using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)

T = Float64

L0_start = 2^3
L0_rate = 2
test_num = 13
point_density = 2^4

tp = TestPolation(; hann_length=3, herm_length=3)

dir1_vec = ["./fft", "./fir", "./loglog", "./asy", "./aaa"]
dir2_vec = ["./type1-schwartz", "./type2-rt_1", "./type3-rt_2", "./type4-rt+s",
            "./type5-logrt", "./type6-logrt+rt+s", "./type7-drt", "./type8-drt+rt+s",
            "./type9-drt+logrt+rt+s", "./type10-drts"]

func_type_vec = [SchwartzFunc(T),
                 RationalFuncPolesRepresent(T; norder=1, npole=[1]),
                 RationalFuncPolesRepresent(T; norder=2, npole=[2, 1]),
                 MixedFunc(T; swf=SchwartzFunc(T),
                           rtfpr=RationalFuncPolesRepresent(T; norder=1, npole=[1])),
                 LogRationalFunc(T; d=1),
                 MixedFunc(T; swf=SchwartzFunc(T),
                           rtfpr=RationalFuncPolesRepresent(T; norder=1, npole=[1]),
                           logrtf=LogRationalFunc(T; d=1)),
                 DRationdlFunc([0.732], T[]),
                 MixedFunc(T; swf=SchwartzFunc(T), drtf=DRationdlFunc([0.732], T[]),
                           rtfpr=RationalFuncPolesRepresent(T; norder=2, npole=[2, 1])),
                 MixedFunc(T; swf=SchwartzFunc(T), drtf=DRationdlFunc([0.732], T[]),
                           logrtf=LogRationalFunc(T; d=1),
                           rtfpr=RationalFuncPolesRepresent(T; norder=2, npole=[2, 1])),
                 DRationdlFunc([0.5, 0.9, 1.7], [0.7, 1.2, 3.8])]

tdm_vec = [TestNoDeMode(), TestNoDeMode(), TestLogLog(; mode_length=10),
           TestAsy(; mode_length=1//10, degree=12, order0=1/3, d=1/3),
           TestAAA(; max_degree=20)]
trans_vec = [FFTTrans(), FIRTrans(), FIRTrans(), FIRTrans(), FIRTrans()]
