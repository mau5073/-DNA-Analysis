import java.util.*;


public class CodonApp
{
    public static void main(String [] args){
        CodonTable condonTable = new CodonTable();
        ArrayList<String> acidSequence = new ArrayList<>();
        Scanner scan = new Scanner(System.in);
        System.out.println("Please enter a codon");
        while(true){
            String codon = scan.nextLine();
            if (codon.equalsIgnoreCase( "quit")) {
                Iterator<String> iterator = acidSequence.iterator();
                        while(iterator.hasNext())
                        {
                            System.out.print(iterator.next());
                        }
                        System.out.println();

                break;
            }
            String acid = condonTable.getAminoAcidSequence( codon );
            if(acid != null){
                acidSequence.add(acid);
            }


            System.out.println("Please enter a codon");
        }
    }
}
