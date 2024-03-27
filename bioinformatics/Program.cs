
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using static System.Formats.Asn1.AsnWriter;

class AlignString 
{ 
    List<String> listOfStrings = new List<String>();
    public AlignString() 
    { 
    
    }
    public void addString(string stringToAdd) 
    { 
        listOfStrings.Add(stringToAdd);
    }
    public int getLength()
    {
        return listOfStrings[0].Length;
    }
    public List<String> getListOfStrings()
    {
        return listOfStrings;
    }



};
class Node
{
    public string Name; 
    public Node Left; // The left child of the node
    public Node Right; // The right child of the node
    public double Distance; // The distance from the parent node
    public AlignString profile;//alignment strings of tree calculate it in every step :(((
    // A constructor to create a leaf node
    public Node(string name ,AlignString profile )
    {
        Name = name;
        Left = null;
        Right = null;
        Distance = 0.0;
        this.profile = profile;
    }

    // A constructor to create an internal node
    public Node(string name,Node left, Node right, double distance, AlignString profile)
    {
        Name = name;
        Left = left;
        Right = right;
        Distance = distance;
        this.profile = profile;

    }


}
public class bioInformatics
{
    public static void Main()
    {

        Dictionary<String, int> BLOSUM62INDEX = new Dictionary<String, int>();
        int[,] BLOSUM62 = new int[24, 24];
        //Console.WriteLine(Directory.GetCurrentDirectory());
        var data = new Dictionary<string, string>();
        int gapValue = -1;
        Console.Write("Bir sayı girin: ");
        string input = Console.ReadLine();

        int number;
        if (Int32.TryParse(input, out number))
        {
            gapValue = number;
        }
        else
        {
            Console.WriteLine("Geçerli bir sayı girin.");
        }

        //input dosyasındaki aminoasitlerin okunup human->ASDAD gibi data dictionarysine yazılması
        using (var reader = new StreamReader("input.txt"))
        {
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.StartsWith(">"))
                {
                    string key = line.Substring(1);  // Remove '>'
                    string value = reader.ReadLine();
                    data[key] = value;
                }
            }
        }
        //BLOSUM MATRİX AND CORRESPONDİNG İNDEX DİCTİONARY CREATİON
        //BLOSUM62INDEX blosum matrixindeki değerlerin indexlerini tutan dictionary
        //BLOSUM62[BLOSUM62INDEX["E"], BLOSUM62INDEX["V"]] bu formatta blosum62 değerlerine erişilebilir
        using (var reader = new StreamReader("Blosum62.txt"))
        {
            int linecount = 0;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                line = line.Trim();
                if (linecount == 0)
                {
                    string[] indexes = Regex.Split(line, @" +");
                    for (int i = 0; i < indexes.Length; i++)
                    {
                        if (indexes[i] == "*")
                        {
                            BLOSUM62INDEX["-"] = i;
                        }
                        else
                        {
                            BLOSUM62INDEX[indexes[i]] = i;
                        }
                    }
                }
                else
                {
                    line = line.Substring(1);
                    line = line.Trim();
                    int[] blosumvalues = Array.ConvertAll(Regex.Split(line, @" +"), int.Parse);
                    for (int i = 0; i < blosumvalues.Length; i++)
                    {
                        BLOSUM62[linecount - 1, i] = blosumvalues[i];
                    }


                }
                linecount++;
            }
            //Console.WriteLine(BLOSUM62[BLOSUM62INDEX["*"], BLOSUM62INDEX["*"]]);
        }
        // Print the dictionary to verify it's correct
        int indexofdata = 0;
        string[] keysofdata = new string[data.Count];
        foreach (var pair in data)
        {
            keysofdata[indexofdata] = pair.Key;
            indexofdata++;
            Console.WriteLine("Key: " + pair.Key + ", Value: " + pair.Value + ",");
        }
        double[,] SimilarityMatrix = CalculateSimilarity(keysofdata.Length);
        double[,] distances = new double[SimilarityMatrix.GetLength(0), SimilarityMatrix.GetLength(0)];
        for (int i = 0; i < SimilarityMatrix.GetLength(0); i++)
        {
            for (int j = 0; j < SimilarityMatrix.GetLength(0); j++)
            {
                if (SimilarityMatrix[i, j] != 0)
                {
                    distances[i, j] = Math.Round(1.0 - SimilarityMatrix[i, j],2); // Convert similarities to distances
                    distances[j, i] = Math.Round(1.0 - SimilarityMatrix[i, j], 2); // Convert similarities to distances
                }
                else
                {
                    distances[i, j] = 0;
                }
            }
        }
        for (int i = 0; i < distances.GetLength(0); i++)
        {
            for (int j = 0; j < distances.GetLength(1); j++)
            {
                Console.Write(distances[i, j].ToString().PadLeft(5));
            }
            Console.WriteLine();
        }
        BuildTree(distances, keysofdata);

        double[,] CalculateSimilarity(int l ) 
        {
            double[,] SimMatrix = new double[l,l];
            for (int i = 0; i < SimMatrix.GetLength(0); i++)
            {
                for (int j = 0; j < SimMatrix.GetLength(0); j++)
                {
                    if (j<i) 
                    {
                        int NumberOfExactMatches = 0;
                        AlignString strinng1 = new AlignString();
                        AlignString strinng2 = new AlignString();
                        strinng1.addString(data[keysofdata[i]]);
                        strinng2.addString(data[keysofdata[j]]);
                        AlignString resultstringlist = MultipleAlign(strinng1, strinng2, gapValue);
                   
                        for (int k = 0; k < resultstringlist.getListOfStrings()[0].Length; k++)
                        {
                            if (resultstringlist.getListOfStrings()[0][k]== resultstringlist.getListOfStrings()[1][k] ) 
                            { 
                                NumberOfExactMatches++;
                            }
                        }
                        SimMatrix[i,j] = Math.Round( (double)NumberOfExactMatches/ resultstringlist.getListOfStrings()[0].Length,2);
                    }
                }

            }
            return SimMatrix;
        }




        Node BuildTree(double[,] distanceMatrix, string[] sequenceNames)
        {
            int n = sequenceNames.Length;
            var nodes = new List<Node>();
            for (int i = 0; i < n; i++)
            {
                AlignString newstring = new AlignString();
                newstring.addString(data[sequenceNames[i]]);
                nodes.Add(new Node(sequenceNames[i], newstring));
            }

            while (nodes.Count > 1)
            {
                int index1 = 0, index2 = 0;
                double min = double.MaxValue;
                for (int i = 0; i < nodes.Count; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        if (distanceMatrix[i, j] < min)
                        {
                            min = distanceMatrix[i, j];
                            index1 = i;
                            index2 = j;
                        }
                    }
                }
                AlignString alignOfTwo = MultipleAlign(nodes[index1].profile, nodes[index2].profile,gapValue);

                //burada multiple alignment yapıp onu üstteki noda eklemek lazımdır
                var newNode = new Node("(" + nodes[index1].Name + "," + nodes[index2].Name + ")", nodes[index1], nodes[index2], min / 2,alignOfTwo);
                //nodes.Insert
                List<Node> oldlist = new List<Node>();
                oldlist.AddRange(nodes);
                nodes.RemoveAt(Math.Max(index1, index2));
                nodes.RemoveAt(Math.Min(index1, index2));
                nodes.Insert(0,newNode);
                distanceMatrix = RecalculateDistances(distanceMatrix, index1, index2, 0,oldlist,nodes);
            }
            Console.WriteLine("MULTİPLE ALİGNMENT:");
            foreach (string s in nodes[0].profile.getListOfStrings())
            {
                Console.WriteLine(s);
            }
            Console.WriteLine("Tree traversal");
            Console.WriteLine(nodes[0].Name);
            return nodes[0];
        }
        void InOrderTraversal(Node node)
        {
            if (node != null)
            { 
                Console.Write(node.Name+" ");
                InOrderTraversal(node.Left);
               
                InOrderTraversal(node.Right);
            }
        }

        double[,] RecalculateDistances(double[,] distanceMatrix, int x, int y, int newNodeIndex,List<Node> oldList,List< Node > newList)
        {
            int n = distanceMatrix.GetLength(0);
            var newMatrix = new double[n - 1, n - 1];
            for (int i = 0; i < newList.Count; i++)
            {
                
                for (int j = 0; j <= i; j++)
                {
                    int oldindexi = -1;
                    int oldindexj = -1;
                    for (int k = 0; k < oldList.Count; k++)
                    {
                        if (oldList[k].Name == newList[i].Name) { oldindexi = k; }
                        if (oldList[k].Name == newList[j].Name) { oldindexj = k; }

                    }

                    if (oldindexi != -1 && oldindexj != -1 && j < i)
                    {
                        newMatrix[i, j] = distanceMatrix[oldindexi, oldindexj];
                        newMatrix[j, i] = distanceMatrix[oldindexi, oldindexj];

                    }
                    if (i == j && oldindexi != -1 && oldindexj != -1)
                    {
                       
                        newMatrix[i, 0] = (distanceMatrix[x, oldindexj] + distanceMatrix[y, oldindexi]) / 2;
                        newMatrix[0, i] = (distanceMatrix[x, oldindexj] + distanceMatrix[y, oldindexi]) / 2;
                    }
                }
            }
           
            return newMatrix;
        }
         

         AlignString MultipleAlign(AlignString str1, AlignString str2, int gapPenalty)
        {
            int m = str1.getLength();

            int n = str2.getLength();

            // Create DP table
            int[,] dp = new int[m + 1, n + 1];
            //bt için 3 ihtimal var sol=1 diagonal=2 yukarı=3;
            int[,] bt = new int[m + 1, n + 1];
            // Initialize first row and column
            List<string> profile1 = str1.getListOfStrings();
            List<string> profile2 = str2.getListOfStrings();
            for (int k = 1; k <= m; k++)
            {
                dp[k, 0] = dp[k - 1, 0] + gapPenalty;
                bt[k, 0] = 3;
            }
            for (int l = 1; l <= n; l++)
            {
                dp[0, l] = dp[0, l - 1] + gapPenalty;
                bt[0, l] = 1;
            }

            // Fill the DP table
            for (int k = 1; k <= m; k++)
            {
                for (int l = 1; l <= n; l++)
                {
                    int score = 0;

                    //calculating the score
                    for (int i1 = 0; i1 < profile1.Count; i1++)
                    {
                        for (int j1 = i1; j1 < profile2.Count; j1++)
                        {
                            score += BLOSUM62[BLOSUM62INDEX[profile1[i1][k - 1].ToString()], BLOSUM62INDEX[profile2[j1][l - 1].ToString()]];
                        }
                    }


                    dp[k, l] = Math.Max(
                        dp[k - 1, l] + gapPenalty, // Delete
                        Math.Max(
                            dp[k, l - 1] + gapPenalty, // Insert
                            dp[k - 1, l - 1] + score // Match/Mismatch
                        )
                    );
                    if ((dp[k - 1, l - 1] + score) >= (dp[k, l - 1] + gapPenalty) && (dp[k - 1, l - 1] + score) >= (dp[k - 1, l] + gapPenalty))//sola ve diyagonale büyük eşittir ile öncelik veriyorum
                    {
                        bt[k, l] = 2;//diagonal geliyorsa
                    }
                    else if ((dp[k, l - 1] + gapPenalty) >= (dp[k - 1, l] + gapPenalty) && (dp[k, l - 1] + gapPenalty) >= (dp[k - 1, l] + gapPenalty))
                    {
                        bt[k, l] = 1;//soldan geliyorsa
                    }
                    else
                    {
                        bt[k, l] = 3;//yukarıdan geliyorsa
                    }
                }
            }

            // Backtracking


            List<StringBuilder> alignedStr3 = new List<StringBuilder>();
            for (int b = 0; b < profile1.Count + profile2.Count; b++)
            {
                alignedStr3.Add(new StringBuilder());
            }
            int i = m, j = n;
            while (i > 0 || j > 0)
            {
                if (i > 0 && j > 0 && bt[i, j] == 2)//diyagonal
                {
                    for (int k = 0; k < profile1.Count; k++)
                    {
                        alignedStr3[k].Insert(0, profile1[k][i - 1]);
                    }
                    for (int k = profile1.Count; k < profile1.Count + profile2.Count; k++)
                    {
                        alignedStr3[k].Insert(0, profile2[k - profile1.Count][j - 1]);
                    }
                    // alignedStr1.Insert(0, str1[i - 1]);
                    //alignedStr2.Insert(0, str2[j - 1]);

                    i--;
                    j--;
                }
                else if (i > 0 && bt[i, j] == 3)//yukarıdansa
                {
                    for (int k = 0; k < profile1.Count; k++)
                    {
                        alignedStr3[k].Insert(0, profile1[k][i - 1]);
                    }
                    for (int k = profile1.Count; k < profile1.Count+ profile2.Count; k++)
                    {
                        alignedStr3[k].Insert(0, '-');
                    }
                    // alignedStr1.Insert(0, str1[i - 1]);
                    // alignedStr2.Insert(0, '-');
                    i--;
                }
                else//soldansa
                {
                    for (int k = 0; k < profile1.Count; k++)
                    {
                        alignedStr3[k].Insert(0, '-');
                    }
                    for (int k = profile1.Count; k < profile1.Count + profile2.Count; k++)
                    {
                        alignedStr3[k].Insert(0, profile2[k-profile1.Count][j - 1]);
                    }
                    //alignedStr1.Insert(0, '-');
                    //alignedStr2.Insert(0, str2[j - 1]);
                    j--;
                }
            }
            AlignString returnprofile = new AlignString();
            for (global::System.Int32 k = 0; k < alignedStr3.Count; k++)
            {
                returnprofile.addString(alignedStr3[k].ToString());
            }
            return returnprofile;
        }

    }
}

//public class Similaritynode
//{
//    public decimal score;
//    public List<String> alignments = new List<String>();
//    public Similaritynode(decimal score)
//    {
//        this.score = score;
//    }
//}



//(int score, string alignedStr1, string alignedStr2) Align(string str1, string str2, int gapPenalty)
//{
//    int m = str1.Length;
//    int score = 0;
//    int n = str2.Length;

//    // Create DP table
//    int[,] dp = new int[m + 1, n + 1];
//    //bt için 3 ihtimal var sol=1 diagonal=2 yukarı=3;
//    int[,] bt = new int[m + 1, n + 1];
//    // Initialize first row and column

//    for (int k = 1; k <= m; k++)
//    {
//        dp[k, 0] = dp[k - 1, 0] + gapPenalty;
//        bt[k, 0] = 3;
//    }
//    for (int l = 1; l <= n; l++)
//    {
//        dp[0, l] = dp[0, l - 1] + gapPenalty;
//        bt[0, l] = 1;
//    }

//    // Fill the DP table
//    for (int k = 1; k <= m; k++)
//    {
//        for (int l = 1; l <= n; l++)
//        {
//            score = BLOSUM62[BLOSUM62INDEX[str1[k - 1].ToString()], BLOSUM62INDEX[str2[l - 1].ToString()]];

//            dp[k, l] = Math.Max(
//                dp[k - 1, l] + gapPenalty, // Delete
//                Math.Max(
//                    dp[k, l - 1] + gapPenalty, // Insert
//                    dp[k - 1, l - 1] + score // Match/Mismatch
//                )
//            );
//            if ((dp[k - 1, l - 1] + score) >= (dp[k, l - 1] + gapPenalty) && (dp[k - 1, l - 1] + score) >= (dp[k - 1, l] + gapPenalty))//sola ve diyagonale büyük eşittir ile öncelik veriyorum
//            {
//                bt[k, l] = 2;//diagonal geliyorsa
//            }
//            else if ((dp[k, l - 1] + gapPenalty) >= (dp[k - 1, l] + gapPenalty) && (dp[k, l - 1] + gapPenalty) >= (dp[k - 1, l] + gapPenalty))
//            {
//                bt[k, l] = 1;//soldan geliyorsa
//            }
//            else
//            {
//                bt[k, l] = 3;//yukarıdan geliyorsa
//            }
//        }
//    }
//    // Backtracking
//    StringBuilder alignedStr1 = new StringBuilder();
//    StringBuilder alignedStr2 = new StringBuilder();
//    int i = m, j = n;
//    while (i > 0 || j > 0)
//    {
//        if (i > 0 && j > 0 && bt[i, j] == 2)//diyagonal
//        {
//            alignedStr1.Insert(0, str1[i - 1]);
//            alignedStr2.Insert(0, str2[j - 1]);
//            i--;
//            j--;
//        }
//        else if (i > 0 && bt[i, j] == 3)//yukarıdansa
//        {
//            alignedStr1.Insert(0, str1[i - 1]);
//            alignedStr2.Insert(0, '-');
//            i--;
//        }
//        else//soldansa
//        {
//            alignedStr1.Insert(0, '-');
//            alignedStr2.Insert(0, str2[j - 1]);
//            j--;
//        }
//    }

//    return (dp[m, n], alignedStr1.ToString(), alignedStr2.ToString());
//}