
remove=feature_
for i in  "$remove"*;do mv "$i" "${i#"$remove"}";done
ls
