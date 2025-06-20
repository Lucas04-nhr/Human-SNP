cd /mnt/raid6/bacphagenetwork/data/12_plink_Full/results

for i in {1..9};
do
  echo "Renaming file result.P${i}.qassoc..."
  mv "result.P${i}.qassoc" "result.P0${i}.qassoc"
  echo "Renaming file result.P${i}.qassoc done."
  echo "Renaming file result.P${i}.assoc.linear..."
  mv "result.P${i}.assoc.linear" "result.P0${i}.assoc.linear"
  echo "Renaming file result.P${i}.assoc.linear done."
done
