
modle="remove2.65"
if modle=="Gx":
    from InvSolver.Seed_algorithm.Solver import Solver
    from InvSolver.Seed_algorithm.Tools import data_tools
    solver=Solver()
    solver.output_Gx_file(r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\output\topo\a1\topo_CK1_4p3p3p_Ball_D3D4low_seed_result",
                          r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\output\topo\a1\Gx"
                          )
elif modle=="remove2.65":
    from InvSolver.Seed_algorithm.Tools import data_tools
    data_tool=data_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr\paper_dafu\gold.den",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr\paper_dafu\ZZG6.msh",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr\paper_dafu\gold.den")
    data_tool.output_res(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr\paper_dafu\gold.den")