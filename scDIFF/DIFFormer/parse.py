
def parser_add_main_args(parser):
    # dataset
    parser.add_argument('--data_dir', type=str, default='../../output/BoneMarrowA_BoneMarrowB/CACNN_output.h5ad')
    parser.add_argument('--device', type=int, default=0,
                        help='which gpu to use if any (default: 0)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--epochs', type=int, default=1000)
    parser.add_argument('--eval_step', type=int,
                        default=1, help='how often to print')
    parser.add_argument('--cpu', action='store_true')
    parser.add_argument('--runs', type=int, default=1,
                        help='number of distinct runs')
    parser.add_argument('--train_prop', type=float, default=.5,
                        help='training label proportion')
    parser.add_argument('--valid_prop', type=float, default=.25,
                        help='validation label proportion')
    parser.add_argument('--metric', type=str, default='acc', choices=['acc', 'rocauc', 'f1'],
                        help='evaluation metric')

    # hyper-parameter for model arch and training

    parser.add_argument('--dropout', type=float, default=0.0)
    parser.add_argument('--lr', type=float, default=0.001)
    parser.add_argument('--weight_decay', type=float, default=0.01)
    parser.add_argument('--batch_size', type=int, default=10000)
    parser.add_argument('--patience', type=int, default=50)

    # hyper-parameter for DIFFormer
    parser.add_argument('--hidden_channels', type=int, default=64)
    parser.add_argument('--num_layers', type=int, default=4,
                    help='number of layers for deep methods')
    parser.add_argument('--num_heads', type=int, default=1)
    parser.add_argument('--alpha', type=float, default=0.5, help='weight for residual link')
    parser.add_argument('--use_residual', action='store_true', help='use residual link for each GNN layer', default=True)
    parser.add_argument('--use_bn', action='store_true', help='use layernorm', default=True)
    parser.add_argument('--use_graph', action='store_true', help='use pos emb', default=True)
    parser.add_argument('--use_weight', action='store_true', help='use weight for GNN convolution', default=True)
    parser.add_argument('--kernel', type=str, default='simple', choices=['simple', 'sigmoid'])


    parser.add_argument("--train_name_list", type=str, nargs='+', default=["BoneMarrow_62216"])
    parser.add_argument("--test_name", type=str, nargs='+', default=["Liver_62016"])
    parser.add_argument("--sample_ratio", type=float, default=0.1)
    parser.add_argument("--edge_ratio", type=float, default=0.0)
    parser.add_argument("--save_path", type=str, default="../../output")
    parser.add_argument("--save_name", type=str, default="BoneMarrowB_liver")

    parser.add_argument("--save_unknown", action='store_true', help='remove unknown cell type', default=False)
    parser.add_argument("--save_rare", action='store_true', help='remove cell types with numbers less than 10', default=False)
    parser.add_argument("--no_smote", action='store_true', help='Up-sampling with smote', default=False)
    parser.add_argument("--get_weight", action='store_true', help='get weight of GraphTransformer', default=False)


