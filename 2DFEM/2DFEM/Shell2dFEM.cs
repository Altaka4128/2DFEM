using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2DFEM
{
    class Shell2dFEM : FEM
    {
        const int ElementNodeNum = 8;
        private int NodeNum = 0;   // 節点数
        public override List<FEM_Element> Elems
        {
            get;
            protected set;
        }

        public Shell2dFEM(FEMData2d data)
        {
            // 節点数を格納する
            NodeNum = data.nodes.Count;

            // 要素の形式を変換して格納する
            Elems = new List<FEM_Element>();
            {
                for (int i = 0; i < data.elems.Count; i++)
                {
                    Node[] nodes = new Node[ElementNodeNum];
                    nodes[0].No = data.elems[i].NodeNo1;
                    nodes[0].Point = data.nodes[nodes[0].No - 1].Point;
                    nodes[1].No = data.elems[i].NodeNo2;
                    nodes[1].Point = data.nodes[nodes[1].No - 1].Point;
                    nodes[2].No = data.elems[i].NodeNo3;
                    nodes[2].Point = data.nodes[nodes[2].No - 1].Point;
                    nodes[3].No = data.elems[i].NodeNo4;
                    nodes[3].Point = data.nodes[nodes[3].No - 1].Point;
                    nodes[4].No = data.elems[i].NodeNo5;
                    nodes[4].Point = data.nodes[nodes[4].No - 1].Point;
                    nodes[5].No = data.elems[i].NodeNo6;
                    nodes[5].Point = data.nodes[nodes[5].No - 1].Point;
                    nodes[6].No = data.elems[i].NodeNo7;
                    nodes[6].Point = data.nodes[nodes[6].No - 1].Point;
                    nodes[7].No = data.elems[i].NodeNo8;
                    nodes[7].Point = data.nodes[nodes[7].No - 1].Point;

                    int materialNo = data.elems[i].MaterialNo - 1;
                    double young = data.materials[materialNo].Young;
                    double poisson = data.materials[materialNo].Poisson;
                    double thickness = data.materials[materialNo].Thickness;

                    Elems.Add(new Shell2dElement(nodes, thickness, young, poisson));
                }
            }

            // 拘束条件の形式を変換して格納する
            // 変位
            List<double> disp = new List<double>();
            List<double> force = new List<double>();
            List<bool> constraint = new List<bool>();
            for (int i = 0; i < data.nodes.Count; i++)
            {
                disp.Add(data.nodes[i].Displacement.X);
                disp.Add(data.nodes[i].Displacement.Y);
                force.Add(data.nodes[i].Force.X);
                force.Add(data.nodes[i].Force.Y);
                constraint.Add(data.nodes[i].ConstraintX);
                constraint.Add(data.nodes[i].ConstraintY);
            }
            DenseVector dispVector = DenseVector.OfArray(disp.ToArray());
            DenseVector forceVector = DenseVector.OfArray(force.ToArray());
            setBoundaryCondition(dispVector, forceVector, constraint);
        }

        // 境界条件を設定する
        public void setBoundaryCondition(DenseVector dispvector, DenseVector forcevector, List<bool> constraint)
        {
            DispVector = dispvector;
            ForceVector = forcevector;
            Constraint = constraint;
        }

        // Kマトリックスを作成する
        private DenseMatrix makeKMatrix()
        {
            // 例外処理
            if (NodeNum <= 0 || Elems == null || Constraint.Count != NodeNum * 2)
            {
                return null;
            }

            DenseMatrix kMatrix = DenseMatrix.Create(NodeNum * 2, NodeNum * 2, 0.0);

            // 各要素のKeマトリックスを計算し、Kマトリックスに統合する
            for (int i = 0; i < Elems.Count; i++)
            {
                Console.WriteLine("要素" + (i + 1).ToString());
                DenseMatrix keMatrix = Elems[i].makeKeMatrix();

                for (int r = 0; r < keMatrix.RowCount; r++)
                {
                    int rt = (Elems[i].Nodes[r / 2].No - 1) * 2 + r % 2;
                    for (int c = 0; c < keMatrix.ColumnCount; c++)
                    {
                        int ct = (Elems[i].Nodes[c / 2].No - 1) * 2 + c % 2;
                        kMatrix[rt, ct] += keMatrix[r, c];
                    }
                }
            }

            Console.WriteLine("Kマトリックス");
            Console.WriteLine(kMatrix);

            // 境界条件を考慮して修正する
            ForceVector = ForceVector - kMatrix * DispVector;
            for (int i = 0; i < Constraint.Count; i++)
            {
                if (Constraint[i] == true)
                {
                    for (int j = 0; j < kMatrix.ColumnCount; j++)
                    {
                        kMatrix[i, j] = 0.0;
                    }
                    for (int k = 0; k < kMatrix.RowCount; k++)
                    {
                        kMatrix[k, i] = 0.0;
                    }
                    kMatrix[i, i] = 1.0;

                    ForceVector[i] = DispVector[i];
                }
            }

            Console.WriteLine("Kマトリックス(境界条件考慮)");
            Console.WriteLine(kMatrix);
            Console.WriteLine("荷重ベクトル(境界条件考慮)");
            Console.WriteLine(ForceVector);

            return kMatrix;
        }

        public override void Analysis()
        {
            DenseMatrix kMatrix = makeKMatrix();

            // 変位を計算する
            DispVector = (DenseVector)(kMatrix.Inverse().Multiply(ForceVector));
            Console.WriteLine("変位ベクトル");
            Console.WriteLine(DispVector);
        }
    }
}
