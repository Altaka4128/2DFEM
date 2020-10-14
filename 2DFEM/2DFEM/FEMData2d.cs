using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace _2DFEM
{
    public struct Material
    {
        public double Young;       // ヤング率
        public double Poisson;     // ポアソン比
        public double Thickness;   // 厚み
    }

    public struct Node
    {
        public int                   No;
        public System.Windows.Vector Point;
        public bool                  ConstraintX;
        public bool                  ConstraintY;
        public System.Windows.Vector Displacement;
        public System.Windows.Vector Force;
    }

    public struct Element
    {
        public int NodeNo1;
        public int NodeNo2;
        public int NodeNo3;
        public int NodeNo4;
        public int NodeNo5;
        public int NodeNo6;
        public int NodeNo7;
        public int NodeNo8;
        public int MaterialNo;
    }

    public struct ResultNode
    {
        public System.Windows.Vector Displacement;
        public System.Windows.Vector Force;
        public DenseVector StressVector;
        public DenseVector PriStressVector;
        public double MisesStress;
    }

    class FEMData2d
    {
        public List<Material> materials = new List<Material>();
        public List<Node> nodes = new List<Node>();
        public List<Element> elems = new List<Element>();
        public List<ResultNode> resultnodes = new List<ResultNode>();

        public void ReadCSVFile1d(String filename)
        {
            /*ファイルフォーマット*/
            //材料
            //Material,Young[Mpa],Poisson,Thickness[mm],,
            //1,210000,0,0.01
            //︙
            //節点
            //Node,X[mm],Y[mm], ConstraintX[mm], ConstraintY[mm], DisplacementX[mm], DisplacementY[mm], ForceX[N], ForceY[N]
            //1,0,0, True, True,0,0,0,0
            //2,5,0, False, False,0,0,0,0
            //︙
            //要素
            //Element,Node1,Node2,Node3,Node4,Material
            //1,1,2,13,12,1
            //2,2,3,14,13,1
            //︙

            var fs = new FileStream(filename, FileMode.Open);
            StreamReader sr = new StreamReader(fs);

            List<String> lines = new List<string>();
            while (!sr.EndOfStream)
            {
                lines.Add(sr.ReadLine());
            }

            // 材料情報の読み込み
            // 材料情報の文字部を読み飛ばす
            int lineNo = 0;
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = true;
                for (int i = 0; i < values.Length; i++)
                {
                    double j = 0;
                    bool result = double.TryParse(values[i], out j);

                    // 文字が含まれている場合、数値まで読み飛ばす
                    if (result == false)
                    {
                        lineNo++;
                        break;
                    }
                    else
                    {
                        headerFlg = false;
                    }
                }
                if (headerFlg == false)
                {
                    lineNo--;
                    break;
                }
            }
            // 材料情報の読み込み
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = false;
                // 最初に文字が含まれている場合、終了する
                double j = 0;
                bool result = double.TryParse(values[0], out j);
                if (result == false && values[0] != "")
                {
                    headerFlg = true;
                    break;
                }
                if (headerFlg == true)
                {
                    break;
                }
                // 材料情報を読み込む
                Material material = new Material();
                material.Young = double.Parse(values[1]);
                material.Poisson = double.Parse(values[2]);
                material.Thickness = double.Parse(values[3]);
                materials.Add(material);
                lineNo++;
            }

            // 節点情報の読み込み
            // 節点情報の文字部を読み飛ばす
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = true;
                for (int i = 0; i < values.Length; i++)
                {
                    double j = 0;
                    bool result = double.TryParse(values[i], out j);

                    // 文字が含まれている場合、数値まで読み飛ばす
                    if (result == false)
                    {
                        lineNo++;
                        break;
                    }
                    else
                    {
                        headerFlg = false;
                    }
                }
                if (headerFlg == false)
                {
                    lineNo--;
                    break;
                }
            }
            // 節点情報の読み込み
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = false;
                // 最初に文字が含まれている場合、終了する
                double j = 0;
                bool result = double.TryParse(values[0], out j);
                if (result == false && values[0] != "")
                {
                    headerFlg = true;
                    break;
                }
                if (headerFlg == true)
                {
                    break;
                }
                // 節点情報を読み込む
                Node node = new Node();
                node.No = int.Parse(values[0]);
                node.Point.X = double.Parse(values[1]);
                node.Point.Y = double.Parse(values[2]);
                node.ConstraintX = Convert.ToBoolean(values[3]);
                node.ConstraintY = Convert.ToBoolean(values[4]);
                node.Displacement.X = double.Parse(values[5]);
                node.Displacement.Y = double.Parse(values[6]);
                node.Force.X = double.Parse(values[7]);
                node.Force.Y = double.Parse(values[8]);
                nodes.Add(node);
                lineNo++;
            }

            // 要素情報の読み込み
            // 要素情報の文字部を読み飛ばす
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = true;
                for (int i = 0; i < values.Length; i++)
                {
                    int j = 0;
                    bool result = int.TryParse(values[i], out j);

                    // 文字が含まれている場合、数値まで読み飛ばす
                    if (result == false)
                    {
                        lineNo++;
                        break;
                    }
                    else
                    {
                        headerFlg = false;
                    }
                }
                if (headerFlg == false)
                {
                    lineNo--;
                    break;
                }
            }
            // 要素情報の読み込み
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = false;
                // 最初に文字が含まれている場合、終了する
                double j = 0;
                bool result = double.TryParse(values[0], out j);
                if (result == false && values[0] != "")
                {
                    headerFlg = true;
                    break;
                }
                if (headerFlg == true)
                {
                    break;
                }
                // 要素情報を読み込む
                Element elem = new Element();
                elem.NodeNo1 = int.Parse(values[1]);
                elem.NodeNo2 = int.Parse(values[2]);
                elem.NodeNo3 = int.Parse(values[3]);
                elem.NodeNo4 = int.Parse(values[4]);
                elem.MaterialNo = int.Parse(values[5]);
                elems.Add(elem);
                lineNo++;
            }

            fs.Close();
        }

        public void ReadCSVFile2d(String filename)
        {
            /*ファイルフォーマット*/
            //材料
            //Material,Young[Mpa],Poisson,Thickness[mm],,
            //1,210000,0,0.01
            //︙
            //節点
            //Node,X[mm],Y[mm], ConstraintX[mm], ConstraintY[mm], DisplacementX[mm], DisplacementY[mm], ForceX[N], ForceY[N]
            //1,0,0, True, True,0,0,0,0
            //2,5,0, False, False,0,0,0,0
            //︙
            //要素
            //Element,Node1,Node2,Node3,Node4,Node5,Node6,Node7,Node8,Material
            //1,1,2,13,12,56,57,58,59,1
            //︙

            var fs = new FileStream(filename, FileMode.Open);
            StreamReader sr = new StreamReader(fs);

            List<String> lines = new List<string>();
            while (!sr.EndOfStream)
            {
                lines.Add(sr.ReadLine());
            }

            // 材料情報の読み込み
            // 材料情報の文字部を読み飛ばす
            int lineNo = 0;
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = true;
                for (int i = 0; i < values.Length; i++)
                {
                    double j = 0;
                    bool result = double.TryParse(values[i], out j);

                    // 文字が含まれている場合、数値まで読み飛ばす
                    if (result == false)
                    {
                        lineNo++;
                        break;
                    }
                    else
                    {
                        headerFlg = false;
                    }
                }
                if (headerFlg == false)
                {
                    lineNo--;
                    break;
                }
            }
            // 材料情報の読み込み
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = false;
                // 最初に文字が含まれている場合、終了する
                double j = 0;
                bool result = double.TryParse(values[0], out j);
                if (result == false && values[0] != "")
                {
                    headerFlg = true;
                    break;
                }
                if (headerFlg == true)
                {
                    break;
                }
                // 材料情報を読み込む
                Material material = new Material();
                material.Young = double.Parse(values[1]);
                material.Poisson = double.Parse(values[2]);
                material.Thickness = double.Parse(values[3]);
                materials.Add(material);
                lineNo++;
            }

            // 節点情報の読み込み
            // 節点情報の文字部を読み飛ばす
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = true;
                for (int i = 0; i < values.Length; i++)
                {
                    double j = 0;
                    bool result = double.TryParse(values[i], out j);

                    // 文字が含まれている場合、数値まで読み飛ばす
                    if (result == false)
                    {
                        lineNo++;
                        break;
                    }
                    else
                    {
                        headerFlg = false;
                    }
                }
                if (headerFlg == false)
                {
                    lineNo--;
                    break;
                }
            }
            // 節点情報の読み込み
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = false;
                // 最初に文字が含まれている場合、終了する
                double j = 0;
                bool result = double.TryParse(values[0], out j);
                if (result == false && values[0] != "")
                {
                    headerFlg = true;
                    break;
                }
                if (headerFlg == true)
                {
                    break;
                }
                // 節点情報を読み込む
                Node node = new Node();
                node.No = int.Parse(values[0]);
                node.Point.X = double.Parse(values[1]);
                node.Point.Y = double.Parse(values[2]);
                node.ConstraintX = Convert.ToBoolean(values[3]);
                node.ConstraintY = Convert.ToBoolean(values[4]);
                node.Displacement.X = double.Parse(values[5]);
                node.Displacement.Y = double.Parse(values[6]);
                node.Force.X = double.Parse(values[7]);
                node.Force.Y = double.Parse(values[8]);
                nodes.Add(node);
                lineNo++;
            }

            // 要素情報の読み込み
            // 要素情報の文字部を読み飛ばす
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = true;
                for (int i = 0; i < values.Length; i++)
                {
                    int j = 0;
                    bool result = int.TryParse(values[i], out j);

                    // 文字が含まれている場合、数値まで読み飛ばす
                    if (result == false)
                    {
                        lineNo++;
                        break;
                    }
                    else
                    {
                        headerFlg = false;
                    }
                }
                if (headerFlg == false)
                {
                    break;
                }
            }
            // 要素情報の読み込み
            while (lines.Count != lineNo)
            {
                String[] values = lines[lineNo].Split(',');
                bool headerFlg = false;
                // 最初に文字が含まれている場合、終了する
                double j = 0;
                bool result = double.TryParse(values[0], out j);
                if (result == false && values[0] != "")
                {
                    headerFlg = true;
                    break;
                }
                if (headerFlg == true)
                {
                    break;
                }
                // 要素情報を読み込む
                Element elem = new Element();
                elem.NodeNo1 = int.Parse(values[1]);
                elem.NodeNo2 = int.Parse(values[2]);
                elem.NodeNo3 = int.Parse(values[3]);
                elem.NodeNo4 = int.Parse(values[4]);
                elem.NodeNo5 = int.Parse(values[5]);
                elem.NodeNo6 = int.Parse(values[6]);
                elem.NodeNo7 = int.Parse(values[7]);
                elem.NodeNo8 = int.Parse(values[8]);
                elem.MaterialNo = int.Parse(values[9]);
                elems.Add(elem);
                lineNo++;
            }

            fs.Close();
        }
        public void WriteCSVFile(String filename)
        {
            // 例外処理
            if(resultnodes == null)
            {
                return;
            }
            
            var fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("Node,DisplacementX,DisplacementY,StressX,StressY,StressXY,PriStress1,PriStress2,MisesStress");

            for (int i = 0; i < resultnodes.Count; i++)
            {
                sw.WriteLine((i + 1) + "," + 
                             resultnodes[i].Displacement.X + "," +
                             resultnodes[i].Displacement.Y + "," +
                             resultnodes[i].StressVector[0] + "," +
                             resultnodes[i].StressVector[1] + "," +
                             resultnodes[i].StressVector[2] + "," +
                             resultnodes[i].PriStressVector[0] + "," +
                             resultnodes[i].PriStressVector[1] + "," +
                             resultnodes[i].MisesStress);
            }

            sw.Close();
            fs.Close();
        }
    }
}
