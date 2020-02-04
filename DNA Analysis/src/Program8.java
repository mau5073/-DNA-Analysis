/**
 * Class Representing Program8.
 * Takes user given file and reads DNA protein encoding regions from letter sequences within file
 * @author Marshall UTISS
 */

import jdk.nashorn.internal.runtime.arrays.ArrayIndex;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;
import java.io.*;

public class Program8  {

// Converts X's and N's within data set to A's
    public static String aConvert( String word ) {

        char[] wordchars = word.toCharArray();
        for (int p = 0; p < word.length(); p++) {
            char   n = word.charAt(p);

            if (n == 'X' || n == 'N') {
                wordchars[p] = 'A';
            }


        }
        word =  String.valueOf(wordchars);

        return word;
    }
  /**  public static String sconvert(String word){

        for(int p = 0; p+4 < word.length(); p++){
            int x = p;
            int y = p+4;
            String n = word.substring(x,y);
          if(n.equals("stop")){
              word.replaceAll(n,"");
            }
        }
        return word;
   */
  public static String mConvert(String word){
        char[] mchars = word.toCharArray();
      for (int p = 0; p < word.length(); p++) {
          char   n = word.charAt(p);
          if (n == 'A') {
              mchars[p] = 'U';
          }
          if (n == 'T') {
              mchars[p] = 'A';
          }
          if (n == 'C') {
              mchars[p] = 'G';
          }
          if (n == 'G') {
              mchars[p] = 'C';
          }

      }
      word =  String.valueOf(mchars);

        return word;
  }

  public static ArrayList<String> aminoAcidSequence(String word) {
      CodonTable condonTable = new CodonTable();
      ArrayList<String> acidSequence = new ArrayList<>();
      for (int i = 0; i+2 < word.length();i++ ) {
          int x = i;
          int y = i + 3;
          String start = word.substring(x,y);
          String acid = condonTable.getAminoAcidSequence( start );

          if(acid.equals("stop")|| acid.equals("TAG") || acid.equals("TAA")){
             acid = acid + System.lineSeparator();
          }
          if (acid != null){
              acidSequence.add(acid);
          }
        i = i+2;
      }
          return acidSequence;

  }



    /**
     * Searches for ATG within letter sequence then looks for one of possible endings to create a
     * protein regions inside an array. The array in this case is a special java array from the
     * ArrayList class. This was done because the number of protein regions that would be found
     * was unknown when making the array and to avoid having null values a dynamic array was needed.
     */
    public static ArrayList<String> atgSearch(String word ) {

        ArrayList<String> regions = new ArrayList<String>();
        String region="";

        int l = 0;
        for (int p = 0; p < word.length() & p+3  < word.length(); p=p+3) {
         boolean   ending= false;

                   int x = p;
                   int y = p + 3;
                   String start = word.substring(x, y);


            if (start.equals("ATG")) {
                int g = y;

                while (g < word.length() & g+3 < word.length() & ending!=true ) {

                    int s = g + 3;
                    int z = g;

                    String end = word.substring(z, s);

                    if (end.equals("TAA") || end.equals("TGA") || end.equals("TAG" )) {
                        region = word.substring(x, s);
                        regions.add(region);
                        l++;
                        ending = true;
                        p = s-3;
                    }
                    if(s+5==word.length()){
                        region = word.substring(x, s+5);
                        regions.add(region+"A"+"TGA");
                        l++;
                        ending = true;
                        p = s-3;
                    }

                    g=g+3;

                }

            }
        }
        return regions;
    }
    // Searches for headers inside file and stores all headers from file into an array
    public static String[] header(String word ) {
        int count = 0;
        char  op = '>';
        for (int i = 0; i < word.length(); i++){
            if(word.charAt(i)==op){
                count ++;
            }
        }
        String headerArray[] = new String[count];
        String header = "";
        boolean end = false;
        int l =0;
        for (int p = 0; p < word.length() && p+3  < word.length(); p++) {
                end = false;
            int x = p;
            String start = word.substring(x, x+1);
                int t = x;
                while(start.equals(">") && t < word.length()&& t+4 < word.length() && end!=true){
                    int o = t;
                    String hEnd = word.substring(o,o+4);

                    if(hEnd.equals("2004") ){
                        end = true;
                        header = word.substring(x,o+4);
                        headerArray[l] = header;
                        l++;
                    }
                        t++;
            }
                p = t;
            }
        return headerArray;
    }
// finds the data between headers  used to search for protein regions stores each into an array

    public static  String[] rData(String word ) {
        int count = 0;
        char  op = '>';
        for (int i = 0; i < word.length(); i++){
            if(word.charAt(i)==op){
                count ++;
            }
        }
        String realData[] = new String[count];
        String data = "";
        boolean end = false;
        int l =0;



        for (int p = 0; p < word.length() && p+4  < word.length(); p++) {
            end = false;
            int x = p+4;
            String start = word.substring(p, x);
            int t = p;
            while(start.equals("2004") && t < word.length()&& t+4 < word.length() && end!=true){
                int o = t;
                String hEnd = word.substring(o,o+1);


                if(hEnd.equals(">")  ){
                    end = true;
                    data = word.substring(p+4,o);
                    realData[l] = data;

                    l++;
                }
                if(o+5 == word.length()){
                    end = true;
                    data = word.substring(p+4,o+5);
                    realData[l] = data;

                    l++;
                }
                t++;
            }
            p = t;
        }
        return realData;
    }

    public static void main(String[] args){

        Scanner keyboard = new Scanner(System.in);
        System.out.print("Please enter the full path to the fasta file.");
        String filename = keyboard.nextLine();
        System.out.print(filename);
        System.out.println("");
        System.out.println("");
// reads user entered data
        File file = new File(filename);
        Scanner inputFile = null;
        try {
            inputFile = new Scanner(file);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        String data = "";

        while (inputFile.hasNext())
        {
             data += inputFile.next();

        }
        inputFile.close();

        System.out.println("What type of data would you like to see?");
        System.out.println("press 1 for PCR's, 2 for mRNA, 3 for Amino Acid Sequences");
        int userchoice = keyboard.nextInt();
        while(userchoice > 3 || userchoice < 1){
            System.out.println( "Please choose a number 1-3");
             userchoice = keyboard.nextInt();
        }
        switch(userchoice) {
            case 1:


// calls header function on data and sets retuned information into headerarray
                String headerArray[] = header(data);
//calls rdata function to limit data to data excluding headers storing them in an array
                String dataArray[] = rData(data);
// prints header data followed by each headers accompanying datas protein regions
                for (int i = 0; i < headerArray.length; i++) {
                    System.out.println(headerArray[i] + "");
                    data = aConvert(dataArray[i]);
                    ArrayList<String> regions = atgSearch(data);
                    for (int b = 0; b < regions.size(); b++) {
                        System.out.println(regions.get(b));
                        System.out.println("");
                    }
                }
            case 2:
                // calls header function on data and sets retuned information into headerarray
                String MheaderArray[] = header(data);
//calls rdata function to limit data to data excluding headers storing them in an array
                String MdataArray[] = rData(data);
                for (int i = 0; i < MheaderArray.length; i++) {
                    System.out.println(MheaderArray[i] + "");
                    data = aConvert(MdataArray[i]);
                    ArrayList<String> regions = atgSearch(data);
                    for(int m = 0; m < regions.size(); m++){
                        String mregion = mConvert(regions.get(m));
                        ArrayList<String> mregions = new ArrayList<>();
                        mregions.add(mregion);


                    for (int b = 0; b < mregions.size(); b++) {
                        System.out.println(mregions.get(b));
                        System.out.println("");
                    }
                    }
                }

            case 3:
                // calls header function on data and sets retuned information into headerarray
                String AheaderArray[] = header(data);
//calls rdata function to limit data to data excluding headers storing them in an array
                String AdataArray[] = rData(data);
              //  StringBuilder amino = new StringBuilder();
                String amino = "";
// prints header data followed by each headers accompanying datas protein regions
                for (int i = 0; i < AheaderArray.length; i++) {
                    System.out.println(AheaderArray[i] + "");
                    data = aConvert(AdataArray[i]);
                    ArrayList<String> regions = atgSearch(data);
                    for( int t = 0; t < regions.size(); t++){
                        ArrayList<String> aminos = aminoAcidSequence(regions.get(t));


                    for (String s : aminos) {

                       // amino.append(s);
                        amino += s;


                    }
                    }

                    System.out.println(amino);
                    amino = "";


                }

        } // leaves user menu

    }
}
