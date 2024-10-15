# SPDX-License-Identifier: BSD-2-Clause

# 7. Simplify climate data for later binning

using JuliaDB
using JuliaDBMeta

@everywhere using Unitful
@everywhere using Unitful.DefaultSymbols
@everywhere using Statistics

cera_era = JuliaDB.load("CERA_ERA_all")

cera_era = @transform cera_era {tmin = min(:t2m, :t2m_1, :t2m_2, :t2m_3, :t2m_4,
                                           :t2m_5, :t2m_6, :t2m_7, :t2m_8,
                                           :t2m_9, :t2m_10, :t2m_11)}
cera_era = @transform cera_era {tmax = max(:t2m, :t2m_1, :t2m_2, :t2m_3, :t2m_4,
                                           :t2m_5, :t2m_6, :t2m_7, :t2m_8,
                                           :t2m_9, :t2m_10, :t2m_11)}
cera_era = @transform cera_era {tmean = mean(uconvert.(K,
                                                       [
                                                           :t2m,
                                                           :t2m_1,
                                                           :t2m_2,
                                                           :t2m_3,
                                                           :t2m_4,
                                                           :t2m_5,
                                                           :t2m_6,
                                                           :t2m_7,
                                                           :t2m_8,
                                                           :t2m_9,
                                                           :t2m_10,
                                                           :t2m_11
                                                       ]))}
cera_era = @transform cera_era {trng = :tmax - :tmin}
cera_era = @transform cera_era {stl1mean = mean(uconvert.(K,
                                                          [
                                                              :stl1,
                                                              :stl1_1,
                                                              :stl1_2,
                                                              :stl1_3,
                                                              :stl1_4,
                                                              :stl1_5,
                                                              :stl1_6,
                                                              :stl1_7,
                                                              :stl1_8,
                                                              :stl1_9,
                                                              :stl1_10,
                                                              :stl1_11
                                                          ]))}
cera_era = @transform cera_era {stl2mean = mean(uconvert.(K,
                                                          [
                                                              :stl2,
                                                              :stl2_1,
                                                              :stl2_2,
                                                              :stl2_3,
                                                              :stl2_4,
                                                              :stl2_5,
                                                              :stl2_6,
                                                              :stl2_7,
                                                              :stl2_8,
                                                              :stl2_9,
                                                              :stl2_10,
                                                              :stl2_11
                                                          ]))}
cera_era = @transform cera_era {stl3mean = mean(uconvert.(K,
                                                          [
                                                              :stl3,
                                                              :stl3_1,
                                                              :stl3_2,
                                                              :stl3_3,
                                                              :stl3_4,
                                                              :stl3_5,
                                                              :stl3_6,
                                                              :stl3_7,
                                                              :stl3_8,
                                                              :stl3_9,
                                                              :stl3_10,
                                                              :stl3_11
                                                          ]))}
cera_era = @transform cera_era {stl4mean = mean(uconvert.(K,
                                                          [
                                                              :stl4,
                                                              :stl4_1,
                                                              :stl4_2,
                                                              :stl4_3,
                                                              :stl4_4,
                                                              :stl4_5,
                                                              :stl4_6,
                                                              :stl4_7,
                                                              :stl4_8,
                                                              :stl4_9,
                                                              :stl4_10,
                                                              :stl4_11
                                                          ]))}
cera_era = @transform cera_era {swvl1mean = mean([
                                                     :swvl1,
                                                     :swvl1_1,
                                                     :swvl1_2,
                                                     :swvl1_3,
                                                     :swvl1_4,
                                                     :swvl1_5,
                                                     :swvl1_6,
                                                     :swvl1_7,
                                                     :swvl1_8,
                                                     :swvl1_9,
                                                     :swvl1_10,
                                                     :swvl1_11
                                                 ])}
cera_era = @transform cera_era {swvl2mean = mean([
                                                     :swvl2,
                                                     :swvl2_1,
                                                     :swvl2_2,
                                                     :swvl2_3,
                                                     :swvl2_4,
                                                     :swvl2_5,
                                                     :swvl2_6,
                                                     :swvl2_7,
                                                     :swvl2_8,
                                                     :swvl2_9,
                                                     :swvl2_10,
                                                     :swvl2_11
                                                 ])}
cera_era = @transform cera_era {swvl3mean = mean([
                                                     :swvl3,
                                                     :swvl3_1,
                                                     :swvl3_2,
                                                     :swvl3_3,
                                                     :swvl3_4,
                                                     :swvl3_5,
                                                     :swvl3_6,
                                                     :swvl3_7,
                                                     :swvl3_8,
                                                     :swvl3_9,
                                                     :swvl3_10,
                                                     :swvl3_11
                                                 ])}
cera_era = @transform cera_era {swvl4mean = mean([
                                                     :swvl4,
                                                     :swvl4_1,
                                                     :swvl4_2,
                                                     :swvl4_3,
                                                     :swvl4_4,
                                                     :swvl4_5,
                                                     :swvl4_6,
                                                     :swvl4_7,
                                                     :swvl4_8,
                                                     :swvl4_9,
                                                     :swvl4_10,
                                                     :swvl4_11
                                                 ])}
cera_era = @transform cera_era {ssrmean = mean([
                                                   :ssr,
                                                   :ssr_1,
                                                   :ssr_2,
                                                   :ssr_3,
                                                   :ssr_4,
                                                   :ssr_5,
                                                   :ssr_6,
                                                   :ssr_7,
                                                   :ssr_8,
                                                   :ssr_9,
                                                   :ssr_10,
                                                   :ssr_11
                                               ])}
cera_era = @transform cera_era {tpmean = mean([
                                                  :tp,
                                                  :tp_1,
                                                  :tp_2,
                                                  :tp_3,
                                                  :tp_4,
                                                  :tp_5,
                                                  :tp_6,
                                                  :tp_7,
                                                  :tp_8,
                                                  :tp_9,
                                                  :tp_10,
                                                  :tp_11
                                              ])}

cera_simple = select(cera_era,
                     (:refval, :date, :tmin, :tmax, :tmean, :trng, :stl1mean,
                      :stl2mean, :stl3mean, :stl4mean, :swvl1mean, :swvl2mean,
                      :swvl3mean, :swvl4mean, :ssrmean, :tpmean))

JuliaDB.save(cera_simple, "CERA_simple")
