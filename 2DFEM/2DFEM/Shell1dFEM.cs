using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2DFEM
{
    class Shell1dFEM : FEM
    {
        public override List<FEM_Element> Elems
        {
            get;
            protected set;
        }
        const int ElementNodeNum = 4;
        private int NodeNum = 0;   // 節点数

        public Shell1dFEM(FEMData2d data)
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

                    int materialNo = data.elems[i].MaterialNo - 1;
                    double young = data.materials[materialNo].Young;
                    double poisson = data.materials[materialNo].Poisson;
                    double thickness = data.materials[materialNo].Thickness;

                    Elems.Add(new Shell1dElement(nodes, thickness, young, poisson));
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

        // 有限要素法を実行する
        public override void Analysis()
        {
            DenseMatrix kMatrix = makeKMatrix();

            // 変位を計算する
            DispVector = (DenseVector)(kMatrix.Inverse().Multiply(ForceVector));
            Console.WriteLine("変位ベクトル");
            Console.WriteLine(DispVector);

            // 各要素の変位、応力ベクトルを計算する
            DenseVector dispElemVector = DenseVector.Create(8, 0.0);
            for (int i = 0; i < Elems.Count; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    dispElemVector[2 * j] = DispVector[2 * (Elems[i].Nodes[j].No - 1)];
                    dispElemVector[2 * j + 1] = DispVector[2 * (Elems[i].Nodes[j].No - 1) + 1];
                }

                Console.WriteLine("要素" + (i + 1).ToString());
                Elems[i].makeStrainVector(dispElemVector);
                Elems[i].makeStressVector();

                // 線形外挿法により各要素の節点応力を求める
                Elems[i].makeNodeStressVector();
            }

            // 各節点の結果を算出する
            List<ResultNode> resultNodes = new List<ResultNode>();
            for (int i = 0; i < NodeNum; i++)
            {
                // 節点応力を計算する
                ResultNode resultNode = new ResultNode();
                resultNode.StressVector = DenseVector.Create(3, 0.0);
                int NodeCount = 0;
                for (int j = 0; j < Elems.Count; j++)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        if (Elems[j].Nodes[k].No == (i + 1))
                        {
                            resultNode.StressVector += Elems[j].NodeStressVector[k];
                            NodeCount++;
                        }
                    }
                }
                resultNode.StressVector /= (double)NodeCount;

                // 変位を格納する
                System.Windows.Vector displacement = new System.Windows.Vector();
                displacement.X = DispVector[2 * i];
                displacement.Y = DispVector[2 * i + 1];
                resultNode.Displacement = displacement;

                // 主応力、ミーゼス応力を計算する
                resultNode.PriStressVector = makePriStress(resultNode.StressVector);
                resultNode.MisesStress = makeMisesStress(resultNode.PriStressVector);

                resultNodes.Add(resultNode);
            }

            // 結果を「output_quad4.csv」で出力する
            FEMData2d data = new FEMData2d();
            data.resultnodes = resultNodes;
            data.WriteCSVFile("output_quad4.csv");

            for (int i = 0; i < NodeNum; i++)
            {
                Console.WriteLine("節点" + (i + 1).ToString());
                Console.WriteLine("節点応力ベクトル");
                Console.WriteLine(resultNodes[i].StressVector);
                Console.WriteLine("節点変位ベクトル");
                Console.WriteLine(resultNodes[i].Displacement + "\n");
                Console.WriteLine("主応力");
                Console.WriteLine(resultNodes[i].PriStressVector);
                Console.WriteLine("ミーゼス応力");
                Console.WriteLine(resultNodes[i].MisesStress + "\n");
            }

        }

        private DenseVector makePriStress(DenseVector stressvector)
        {
            // 例外処理
            if (stressvector == null)
            {
                return null;
            }

            // 主応力を計算する
            double[,] stressTensorArray = new double[2, 2];
            stressTensorArray[0, 0] = stressvector[0];
            stressTensorArray[1, 1] = stressvector[1];
            stressTensorArray[0, 1] = stressvector[2];
            stressTensorArray[1, 0] = stressTensorArray[0, 1];
            DenseMatrix stressTensor = DenseMatrix.OfArray(stressTensorArray);
            var evd = stressTensor.Evd();
            var evdValue = evd.EigenValues;   // 固有値

            Console.WriteLine("主応力");
            Console.WriteLine(evdValue);

            DenseVector evdVector = DenseVector.Create(2, 0.0);
            evdVector[0] = evdValue[0].Real;
            evdVector[1] = evdValue[1].Real;

            return evdVector;

        }

        private double makeMisesStress(DenseVector pristressvector)
        {
            // 例外処理
            if (pristressvector == null)
            {
                return 0.0;
            }

            double misesStress = new double();
            misesStress = (pristressvector[0] * pristressvector[0] + pristressvector[1] * pristressvector[1] +
                           (pristressvector[1] - pristressvector[0]) * (pristressvector[1] - pristressvector[0])) / 2.0;
            misesStress = Math.Sqrt(misesStress);

            Console.WriteLine("ミーゼス応力");
            Console.WriteLine(misesStress + "\n");

            return misesStress;
        }
    }
}
