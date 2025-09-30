# Graph-Regularizer-Survival

This repo contains the official implementation for the paper **[Survival Analysis with Graph-Based Regularization for Predictors](https://link.springer.com/article/10.1007/s12561-025-09483-8)** \
by [Liyan Xie](https://liyanxie.github.io/), Xi He, [Pinar Keskinocak](https://sites.gatech.edu/pinar-keskinocak/) and [Yao Xie](https://www2.isye.gatech.edu/~yxie77/).


We implement a graph-regularized Cox proportional hazards model for variable selection in survival analysis, leveraging prior knowledge of correlations among variables encoded as a graph.


Main code: graph-regularizer-main.R

Supporting functions are provided in: graph-lasso-helper.R


## References

If you find the code useful for your research, please consider citing

```bibtex
@article{XHKX2025survival,
  title={Survival Analysis with Graph-based Regularization for Predictors},
  author={Liyan Xie and Xi He and Pinar Keskinocak and Yao Xie},
  journal={Statistics in Biosciences},
  pages={1--36},
  year={2025},
  publisher={Springer}
}

```

This work is built upon some previous papers which might also interest you:

- Guan Yu and Yufeng Liu. "Sparse regression incorporating graphical structure among predictors." Journal of the American Statistical Association 111.514 (2016): 707-720.
- Jianqing Fan and Runze Li. "Variable selection for Cox's proportional hazards model and frailty model." The Annals of Statistics 30.1 (2002): 74-99.
  
