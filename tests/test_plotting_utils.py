from test_pipeline import create_fake_data
from kinase_msm.series_setup import setup_series_analysis
from kinase_msm.plotting_utils import *
from kinase_msm.fit_transform_kinase_series import fit_pipeline
from kinase_msm.mdl_analysis import Protein,ProteinSeries
from mdtraj.utils.contextmanagers import enter_temp_directory
import os

def test_plotting_utils():
    with enter_temp_directory():
        base_dir = os.path.abspath(os.path.curdir)
        mdl_dir = os.path.join(base_dir,"mdl_dir")
        feature_dir = "feature_dir"
        series_name = "fake_series"
        protein_list = ["kinase_1", "kinase_2"]
        project_dict = {"kinase_1": ["fake_proj1",],
                            "kinase_2": ["fake_proj2"]}
        mdl_params = {'tica__n_components': 1, 'tica__lag_time': 1,
                      'tica__weighted_transform': True, 'tica__shrinkage': 0.01,
                      'cluster__n_clusters': 2,'msm__lag_time': 1,
                      'bootrap__n_samples':1
                      }


        create_fake_data(base_dir, protein_list, project_dict)
        setup_series_analysis(base_dir, mdl_dir, feature_dir,
                                      series_name, protein_list,
                                      project_dict, mdl_params)

        fit_pipeline(base_dir)
        prj = ProteinSeries(os.path.join(mdl_dir,"project.yaml"))

        prt1 = Protein(prj, "kinase_1")
        prt2 = Protein(prj, "kinase_2")

        n_bins = 100

        lin_spaced_tic_dict = global_tic_boundaries([prt1, prt2],
                                                    range(prt1.n_tics_), n_bins)

        def test_bounds():
            locally_calc={}
            for i in range(prt1.n_tics_):
                locally_calc[i] =[]
                global_min = min(min([min(i) for i in prt1.tica_data.values()]),
                    min([min(i) for i in prt2.tica_data.values()]))

                locally_calc[i].append(global_min)

                global_max = max(max([max(i) for i in prt1.tica_data.values()]),
                    max([max(i) for i in prt2.tica_data.values()]))

                locally_calc[i].append(global_max)

            for i in range(prt1.n_tics_):
                assert(lin_spaced_tic_dict[i][0]==locally_calc[i][0])
                assert(lin_spaced_tic_dict[i][-1]==locally_calc[i][-1])
                assert(len(lin_spaced_tic_dict[i])==n_bins)

            return True

        def test_histogram_data():
            H_dict, H_calc, _ = tica_histogram(prj, prt1, [0],
                                               x_array=lin_spaced_tic_dict[0],
                                               n_bins=None)
            assert(len(H_dict.keys()) == prt1.n_states_)
            assert(len(H_calc) == len(lin_spaced_tic_dict[0])-1)
            rnd_state = np.random.randint(0,prt1.n_states_)
            assert(np.allclose(H_dict[rnd_state], np.histogram(prt1.tic_dict[0][rnd_state],
                                                       bins = lin_spaced_tic_dict[0],
                                                       normed=True)[0]))
            return True


        def test_one_dim_free_energy():
            df = one_dim_tic_free_energy(prj, prt1, 0, n_bins=None ,
                        lin_spaced_tic=lin_spaced_tic_dict[0], errorbars=False)

            assert((df.protein_name==prt1.name).all())
            assert((df.mdl_index=="mle").all())

            return True

        assert(test_bounds())
        assert(test_histogram_data())
        assert(test_one_dim_free_energy())


        return