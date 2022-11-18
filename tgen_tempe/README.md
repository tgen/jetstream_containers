# Tgen-Pipelines

Open Container Initiative containers for use at TGen.

## Use Case: 
- Create Containers with Only Providing a Dockerfile
- Images will be created and pushed to:
  - **Artifactory:** [Artifactory Image Registry ](http://pbc-art-prd01:8082/ui/repos/tree/General/)
  - **Github Packages:**  [Github Image Registry ](https://github.com/orgs/tgen/packages?repo_name=jetstream_containers) 
---
## Pipeline Instructions
---
### How are the Container Images Tagged:
    Will use folder path to name and tag containers
    Example:

    Individual_Images/p/python/2.7/Dockerfile => <IMAGE-REPO-URL>/python:2.7
    
    Individual_Images/t/my-terraform-image/12.01/Dockerfile => <IMAGE-REPO-URL>/my-terraform-image:12.01


- When creating a new container image, make sure your folder path follows the structure as shown above to properly tag your image.

**Disclaimer** -> DO NOT CAPITALIZE ANY LETTER WHEN CREATING THE FOLDER PATH. Container images can not have capitalized letters in their name.

---
###  **Create a New Container:**
- 1a.) First You'll Need to Git Clone the Repo:

    ```
    git clone https://gitlab01.tgen.org/oci-containers/external-images/tgen-pipelines.git 
    ```
OR
- 1b.) If you already have the Gitlab repo on your machine. Update your `main` branch before making changes by using the following command:
    ```
    git pull origin main
    ```
---
- 2.) Create New Branch:
    ```
    # To create a new branch use this command:
    git checkout -b dev-agomez
    ```
---
- 3a). Create a Folder
    - Create New Folders for Image Name and Version
    ```
    # I want to create an image called python-flask with version 3.9. The python-flask folder does not exist:

    mkdir -p Individual_Images/p/python-flask/3.9
    touch Individual_Images/p/python-flask/3.9/Dockerfile
    ```
OR
- 3b.) Use Existing Folder but Different Version
    - Use Same Folder name but different version. The python folder already exisit but I'm creating a new version with a different tag
    ```
    # I want to create an image called python with version 3.10
    touch Individual_Images/p/python/3.10/Dockerfile
    ```
---
- 3.) Write out and configure the Dockerfile
---
- 4.) Start Git Workflow for Gitlab
  ```
  # 1.) git add <path-to-your-Dockerfile> 
  git add Individual_Images/p/python/3.10/Dockerfile
  
  # 2.) git commit -m "Description-of-change"
  git commit -m "Added image python 3.10 for dev testing"

  # 3.) git push origin <your-newly-created-branch>
  git push origin dev-agomez
  ```
---
- 5.) Your `git push` command will now trigger the Jenkins Pipeline. This will:
    - Create your Container Image
    - Push it to `Artifactory` -> [Artifactory Image Registry ](http://pbc-art-prd01:8082/ui/repos/tree/General/)

    You can see the status of the build pipeline at:
    - [Jenkins Pipeline Status](https://jenkins01:8443/job/OCI-Containers/job/External-OCI-Pipeline/)
    - [Gitlab Pipeline Status](https://gitlab01.tgen.org/oci-containers/external-images/tgen-pipelines/-/pipelines)
---
- 6.) Test your image and validate until it's functioning properly.
    Remember to update your Gitlab Branch by following the Git workflow for every new change. 
    **Refer to Step #4 for more information**
---
- 7.) Once your image is working properly submit a MR to `main` branch to merge the changes from your branch to main(aka prod branch)
    - On the Gitlab Repo main page:
      - Click to the `"Merge Request"` on the left tab
    - On the `"Merge requests"` page:
      - Click on `"New merge request"` on the top right
      - Click on `"Source Branch"` `"Select source branch"` and select your branch name
      - Leave `"Target Branch"` at `"main"` branch
      - Click `"Compare branches and continue"`
    - On the `"New merge request"` page:
      - `"Title"` = "username-image_added"
      - `"Description"` = Gist for the reasoning of the merge request
      - `"Assignee"` = The user submitting the MR
      - `"Reviewer"` = "Andy Gomez @agomez"
      - For `"Merge options"` check `"Delete source branch when merge request is accepted."`
      - Leave everything else at there Default value
      - Now wait until the MR is approved
---
- 8.) Once Approved, the newly added image will:
    - Be added to the `"main"` branch and it's now protected. **Protected** means that no one can make changes to your image unless they submit a MR and it gets approved.
  - Image will be pushed to Github Packages: [Github Image Registry ](https://github.com/orgs/tgen/packages?repo_name=jetstream_containers)
---
---
### **Update an Existing Container:**
- 1a.) First You'll Need to Git Clone the Repo:
    ```
    git clone https://gitlab01.tgen.org/oci-containers/external-images/tgen-pipelines.git 
    ```

- 1b.) If you already have the Gitlab repo on your machine. Make you update your `main` branch before making changes by using the following command:
    ```
    git pull origin main
    ```
---
- 2.) Create New Branch:
    ```
    # To create a new branch use this command:
    git checkout -b update-python-3.10
    ```
---
- 3.) Go to an existing Dockerfile and make changes and/or update the file 
---
- 4.) Start Git Workflow for Gitlab
  ```
  # 1.) git add <path-to-your-Dockerfile> 
  git add Individual_Images/p/python/3.10/Dockerfile
  
  # 2.) git commit -m "Description-of-chang"
  git commit -m "Updated image python 3.10 for newly added dependencies"

  # 3.) git push origin <your-newly-created-branch>
  git push origin update-python-3.10
  ```
---
- 5.) Once pushed to Gitlab, submit a MR for main branch to apply the updates the image:
    - On the Gitlab Repo main page:
      - Click to the `"Merge Request"` on the left tab
    - On the `"Merge requests"` page:
      - Click on `"New merge request"` on the top right
      - Click on `"Source Branch"` `"Select source branch"` and select your branch name
      - Leave `"Target Branch"` at `"main"` branch
      - Click `"Compare branches and continue"`
    - On the `"New merge request"` page:
      - `"Title"` = "username-image_added"
      - `"Description"` = Gist for the reasoning of the merge request
      - `"Assignee"` = The user submitting the MR
      - `"Reviewer"` = "Andy Gomez @agomez"
      - For `"Merge options"` check `"Delete source branch when merge request is accepted."`
      - Leave everything else at there Default value
      - Now wait until the MR is approved
---
- 6.) Once Approved, the updated image will:
    - Apply new updates to the image on the "main" branch.
    - Will push out updated image to:
      - **Artifactory:** [Artifactory Image Registry ](http://pbc-art-prd01:8082/ui/repos/tree/General/)
      - **Github Packages:**  [Github Image Registry ](https://github.com/orgs/tgen/packages?repo_name=jetstream_containers) 
---
**Disclaimer**:
Update existing images will override any image on both Artfactory and Github Packages. Best Practice to update existing images will be to create a new folder at the version level with "v2" at the end of the tag. For example:

- I want to update python 3.10
    ```
    ls Individual_Images/p/python/

    --Output -> 3.10
    ```
- Create a folder with the same version name but with "-v2" at the end.
    ```
    mkdir Individual_Images/p/python/3.10-v2
    ls Individual_Images/p/python

    --Output -> 3.10 3.10-v2
    ```

- Create Dockerfile for you new version
    ```
    touch Individual_Images/p/python/3.10-v2/Dockerfile
    ```

- Copy over the Dockerfile from Individual_Images/p/python/3.10 into your new v2 folder and follow "Create a New Container" steps from above

---
---
## Local Machine Instructions:
---
### **How to create an image from a Dockerfile:**

    1.) sudo podman login --tls-verify=false pbc-art-prd01.ad.tgen.org

    2.) sudo buildah bud -f {dockerfile)} -t {image name} .

    3.) sudo podman tag {image:tag} pbc-art-prd01.ad.tgen.org/hpc-local-containers/{image name}:{tag}

    4.) sudo podman push pbc-art-prd01.ad.tgen.org/hpc-local-containers/{image name}:{tag}

### **How to pull an image from Artifactory:**
    1.) sudo podman login --tls-verify=false pbc-art-prd01.ad.tgen.org

    2.) sudo podman pull pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/{image name}:{tag}
---
---
## Documentation
- [Jira - Project Pages](https://it-tgen.atlassian.net/jira/software/c/projects/OP/pages)
---
## Org Resources
- [Gitlab Project - OCI Containers ](https://gitlab01.tgen.org/oci-containers/external-images/tgen-pipelines)
- [Artifactory - Private Image Registry](http://pbc-art-prd01:8082/ui/repos/tree/General/)
- [Jenkins - Pipeline](https://jenkins01:8443/job/OCI-Containers/job/External-OCI-Pipeline/)
- [Github Project - OCI Containers](https://github.com/tgen/jetstream_containers/tree/main/tgen_maricopa)
- [Github Packages - Public Image Registry](https://github.com/orgs/tgen/packages?repo_name=jetstream_containers)
---