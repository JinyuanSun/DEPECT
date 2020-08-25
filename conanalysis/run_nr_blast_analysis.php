
<?php

$seqfilename = "uploadseq/seq" . date("Y-m-d-h-i-sa") . ".fasta";

$ic = $_GET['identity'];
$e = $_GET['evalue'];
$score=$_GET['score'];
file_put_contents($seqfilename, $_GET['sequence']."\n");

$command = "python3 depect_ca_server.py -s ".$seqfilename." -db  /ndata/database/blastdb/nr -e ".$e." -ic ".$ic." -nt 4 -sc ".$score." -m blast";
echo exec($command);
sleep(40);
function downtemplateAction($file_name){
    header("Content-type:text/html;charset=utf-8");

    $file_name = iconv("utf-8","gb2312",$file_name);
    $file_sub_path = "/var/www/html/depect/";
    $file_path=$file_sub_path.$file_name;
    if(!file_exists($file_path))
    {
        echo "计算还未完成，请稍后访问 www.lab.cf/depect/".$file_name." 下载！";
        exit;
    }
 
    $fp=fopen($file_path,"r");
    $file_size=filesize($file_path);

    Header("Content-type: application/octet-stream");
    Header("Accept-Ranges: bytes");
    Header("Accept-Length:".$file_size);
    Header("Content-Disposition: attachment; filename=".$file_name);
    $buffer=1024;
    $file_count=0;
    while(!feof($fp) && $file_count<$file_size)
    {
        $file_con=fread($fp,$buffer);
        $file_count+=$buffer;
        echo $file_con;
    }
    fclose($fp);    
}

ob_clean();
ob_end_flush();
downtemplateAction($seqfilename.".ca");

?>
